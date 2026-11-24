from __future__ import annotations

import hashlib
import json
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam

from .cohort_vcf import expand_chroms
from .personal_fasta import reference_fingerprint


ASSEMBLY_QC_SCHEMA = "mfm_personal_fasta_sample_qc_v1"
ASSEMBLY_SOURCE = "assembly_direct"
ALLOWED_BASES = frozenset("ACGTN")


@dataclass
class AssemblyDirectConfig:
    cohort: str
    sample: str
    output_root: Path
    reference_fasta: Path
    hap1_fasta: Path
    hap1_chrom_alias: Path
    hap2_fasta: Path
    hap2_chrom_alias: Path
    chroms: list[str]
    hap1_chain: Path | None = None
    hap2_chain: Path | None = None
    phase_label_status: str = "unverified"
    overwrite: bool = False


def run_assembly_direct(config: AssemblyDirectConfig) -> dict[str, Any]:
    """표준 chromosome-level assembly를 기존 per-chromosome FASTA 계약으로 변환한다.

    HPRC chromAlias에서 요청 chromosome과 정확히 일치하는 단일 contig만 허용한다.
    random/unplaced contig를 이어 붙이는 것은 좌표를 새로 발명하는 작업이므로 이
    adapter에서는 의도적으로 지원하지 않는다.
    """

    started = time.time()
    chroms = expand_chroms(config.chroms)
    if not chroms:
        raise ValueError("assembly-direct chromosome list is empty")
    if config.phase_label_status not in {"unverified", "validated"}:
        raise ValueError("phase_label_status must be unverified|validated")
    _require_fasta(config.reference_fasta, "reference_fasta")

    sample_root = config.output_root.resolve() / config.sample
    individual_root = sample_root / "individual_fasta"
    qc_path = individual_root / "personal_fasta_qc.json"
    source_specs = {
        "HAP1": (config.hap1_fasta, config.hap1_chrom_alias, config.hap1_chain),
        "HAP2": (config.hap2_fasta, config.hap2_chrom_alias, config.hap2_chain),
    }
    source_identity = {
        hap: _source_identity(config.sample, hap, fasta, alias, chain)
        for hap, (fasta, alias, chain) in source_specs.items()
    }
    if qc_path.exists() and not config.overwrite:
        previous = json.loads(qc_path.read_text(encoding="utf-8"))
        if _reusable(previous, config, chroms, source_identity, individual_root):
            return {**previous, "run_status": "SKIP", "message": "existing compatible assembly-direct FASTA"}
        raise FileExistsError(
            f"기존 individual_fasta가 assembly-direct 입력 계약과 다릅니다. "
            f"검토 후 --overwrite를 명시하세요: {individual_root}"
        )

    individual_root.mkdir(parents=True, exist_ok=True)
    hap_results: dict[str, list[dict[str, Any]]] = {}
    hap_summaries: dict[str, dict[str, Any]] = {}
    for hap, (fasta_path, alias_path, _chain_path) in source_specs.items():
        alias = _load_exact_chrom_alias(alias_path, config.sample, hap, chroms)
        rows, summary = _extract_haplotype(
            sample=config.sample,
            hap=hap,
            source_fasta=fasta_path,
            aliases=alias,
            chroms=chroms,
            output_dir=individual_root / hap,
            overwrite=config.overwrite,
        )
        hap_results[hap] = rows
        hap_summaries[hap] = summary

    qc: dict[str, Any] = {
        "schema": ASSEMBLY_QC_SCHEMA,
        "status": "pass",
        "sample": config.sample,
        "cohort": config.cohort,
        "sequence_source": ASSEMBLY_SOURCE,
        "sequence_variant_policy": "not_applicable_assembly",
        "coordinate_space": "assembly_chromosome_contig_renamed_to_chrN",
        "chroms": chroms,
        "reference_fasta": str(config.reference_fasta.resolve()),
        "reference_role": "GRCh38 coordinate/provenance link; not assembly sequence source",
        "reference_fingerprint": reference_fingerprint(config.reference_fasta),
        "chain_policy": "official multi-chain retained as provenance only; not parsed as bcftools single chain",
        "blacklist_policy": "GRCh38 blacklist disabled because official multi-chain projection is not implemented",
        "fragment_policy": "fail_if_exactly_one_chromosome_level_contig_is_not_available",
        "phase_label_status": config.phase_label_status,
        "safe_hap_read_mode": "all_to_both" if config.phase_label_status != "validated" else "validated_labels_allow_opt_in_routing",
        "source_identity": source_identity,
        "hap_results": hap_results,
        "hap_summaries": hap_summaries,
        "individual_fasta_dir": str(individual_root),
        "created_at_unix": time.time(),
        "elapsed_sec": round(time.time() - started, 3),
    }
    _write_json_atomic(qc_path, qc)
    return {**qc, "run_status": "OK", "message": "created"}


def _source_identity(
    sample: str,
    hap: str,
    fasta: Path,
    alias: Path,
    chain: Path | None,
) -> dict[str, Any]:
    _require_fasta(fasta, f"{hap} source FASTA")
    _require_file(alias, f"{hap} chromAlias")
    if chain is not None:
        _require_file(chain, f"{hap} official chain")
    return {
        "sample": sample,
        "haplotype": hap,
        "fasta": _file_identity(fasta),
        "fasta_index": _file_identity(Path(str(fasta) + ".fai")),
        "chrom_alias": _file_identity(alias),
        "official_chain": _file_identity(chain) if chain is not None else None,
        "official_chain_direction": "assembly(target)->GRCh38(query)" if chain is not None else "not_provided",
    }


def _load_exact_chrom_alias(
    path: Path,
    sample: str,
    hap: str,
    chroms: list[str],
) -> dict[str, str]:
    wanted = set(chroms)
    matches: dict[str, list[str]] = {chrom: [] for chrom in chroms}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"invalid chromAlias row: {path}: {line.rstrip()}")
            assembly_name, ucsc_name = fields[0], fields[1]
            if ucsc_name not in wanted:
                continue
            if not assembly_name.startswith(f"{sample}#"):
                raise ValueError(
                    f"chromAlias sample mismatch: expected={sample} hap={hap} assembly={assembly_name}"
                )
            matches[ucsc_name].append(assembly_name)
    missing = [chrom for chrom in chroms if not matches[chrom]]
    duplicate = {chrom: names for chrom, names in matches.items() if len(names) != 1}
    if missing or duplicate:
        raise ValueError(
            "assembly-direct requires exactly one chromosome-level contig per requested chromosome; "
            f"hap={hap} missing={missing} duplicate={duplicate}. random/unplaced contigs are not concatenated"
        )
    return {chrom: matches[chrom][0] for chrom in chroms}


def _extract_haplotype(
    *,
    sample: str,
    hap: str,
    source_fasta: Path,
    aliases: dict[str, str],
    chroms: list[str],
    output_dir: Path,
    overwrite: bool,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    selected_names = set(aliases.values())
    with pysam.FastaFile(str(source_fasta)) as source:
        source_lengths = dict(zip(source.references, source.lengths, strict=True))
        resolved_aliases: dict[str, str] = {}
        for chrom, assembly_name in aliases.items():
            candidates = (assembly_name, assembly_name.rsplit("#", 1)[-1])
            found = [name for name in dict.fromkeys(candidates) if name in source_lengths]
            if len(found) != 1:
                raise ValueError(
                    f"chromAlias contig를 source FASTA에서 유일하게 찾을 수 없습니다: "
                    f"hap={hap} chrom={chrom} candidates={candidates} found={found}"
                )
            resolved_aliases[chrom] = found[0]
        selected_names = set(resolved_aliases.values())
        omitted_names = [name for name in source.references if name not in selected_names]
        omitted_bases = sum(source_lengths[name] for name in omitted_names)
        for chrom in chroms:
            source_name = resolved_aliases[chrom]
            sequence = source.fetch(source_name).upper()
            invalid: dict[str, int] = {}
            for base in sequence:
                if base not in ALLOWED_BASES:
                    invalid[base] = invalid.get(base, 0) + 1
            if invalid:
                raise ValueError(f"assembly FASTA에 A/C/G/T/N 외 base가 있습니다: {hap}/{chrom} {invalid}")
            fasta_path = output_dir / f"{chrom}.fa"
            index_path = Path(str(fasta_path) + ".fai")
            if (fasta_path.exists() or index_path.exists()) and not overwrite:
                raise FileExistsError(f"기존 assembly-direct contig 출력이 있습니다: {fasta_path}")
            partial = output_dir / f".{chrom}.fa.partial.{os.getpid()}"
            with partial.open("w", encoding="ascii") as handle:
                handle.write(f">{chrom}\n")
                for start in range(0, len(sequence), 80):
                    handle.write(sequence[start : start + 80] + "\n")
            os.replace(partial, fasta_path)
            index_path.unlink(missing_ok=True)
            pysam.faidx(str(fasta_path))
            with pysam.FastaFile(str(fasta_path)) as check:
                if tuple(check.references) != (chrom,) or check.get_reference_length(chrom) != len(sequence):
                    raise RuntimeError(f"assembly-direct FASTA index 검증 실패: {fasta_path}")
            rows.append(
                {
                    "sample": sample,
                    "haplotype": hap,
                    "contig": chrom,
                    "source_contig": source_name,
                    "source_length": source_lengths[source_name],
                    "consensus_length": len(sequence),
                    "n_count": sequence.count("N"),
                    "n_fraction": sequence.count("N") / len(sequence) if sequence else 0.0,
                    "invalid_bases": invalid,
                    "fasta": str(fasta_path),
                    "fasta_md5": _md5(fasta_path),
                    "chain": "",
                    "chain_qc": {"status": "not_applicable_assembly_multichain"},
                    "status": "pass",
                }
            )
    return rows, {
        "selected_contigs": len(selected_names),
        "selected_bases": sum(row["consensus_length"] for row in rows),
        "omitted_contigs": len(omitted_names),
        "omitted_bases": omitted_bases,
        "omitted_policy": "random_unplaced_and_nonrequested_contigs_excluded",
    }


def _reusable(
    previous: dict[str, Any],
    config: AssemblyDirectConfig,
    chroms: list[str],
    source_identity: dict[str, Any],
    individual_root: Path,
) -> bool:
    if not (
        previous.get("schema") == ASSEMBLY_QC_SCHEMA
        and previous.get("status") == "pass"
        and previous.get("sample") == config.sample
        and previous.get("sequence_source") == ASSEMBLY_SOURCE
        and previous.get("chroms") == chroms
        and previous.get("phase_label_status") == config.phase_label_status
        and previous.get("source_identity") == source_identity
    ):
        return False
    return all(
        (individual_root / hap / f"{chrom}.fa").is_file()
        and Path(str(individual_root / hap / f"{chrom}.fa") + ".fai").is_file()
        for hap in ("HAP1", "HAP2")
        for chrom in chroms
    )


def _require_fasta(path: Path, label: str) -> None:
    _require_file(path, label)
    _require_file(Path(str(path) + ".fai"), f"{label} index")


def _require_file(path: Path, label: str) -> None:
    if not path.exists() or not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"{label}가 없거나 비었습니다: {path}")


def _file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    return {"path": str(resolved), "size": resolved.stat().st_size, "md5": _md5(resolved)}


def _md5(path: Path, block_size: int = 8 << 20) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f".{path.name}.partial.{os.getpid()}")
    partial.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(partial, path)

