from __future__ import annotations

import csv
import hashlib
import json
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import pysam

from pre_training.preprocess_raw_readset.coordinates import (
    ChainAlignment,
    interval_mapping_coverage,
    load_chain_alignment,
)
from pre_training.preprocess_raw_readset.windowing import deterministic_offset, trim_terminal_n, sequence_annotations
from pre_training.preprocess_raw_readset.window_filter_stats import (
    DROP_REASON_PRIORITY,
    N_FRACTION_DENOMINATOR,
    count_index_exclusions,
    exclusion_reason,
)

from .cohort_vcf import expand_chroms, load_sample_ids
from .sample_vcf import normalize_sequence_variant_policy
from .window_contract import (
    N_VALUE_FIELDS, ensure_window_contract, ensure_window_n_annotations, n_value_metadata, read_windows,
)


DEFAULT_BLACKLIST_BED = "references/hg38_autosomes_blacklist_union.bed"
WINDOW_INDEX_SCHEMA = "mfm_integrated_window_index_v2"
PROCESS_CODE_WINDOW_SCHEMA = "mfm_window_index_v6"
ACGT = frozenset("ACGT")


@dataclass
class WindowIndexConfig:
    """Configuration for process_code-style sample window index generation.

    The blacklist/N/QC policy follows process_code/05_windows.py.  The window
    start policy intentionally follows the integrated preprocessing decision:
    deterministic random offsets by sample/haplotype/contig unless offset_unit=1
    is explicitly requested.
    """

    cohort: str
    output_root: Path
    chroms: list[str]
    sequence_variant_policy: str
    sequence_source: str = "vcf_consensus"
    sample_list: Path | None = None
    samples: list[str] | None = None
    window_size: int = 1 << 20
    offset_unit: int | None = None
    seed: int = 42
    drop_n: bool = False
    blacklist_bed: Path | None = Path(DEFAULT_BLACKLIST_BED)
    disable_reference_blacklist: bool = False
    max_blacklist_frac: float = 1.0
    jobs: int = 8
    overwrite: bool = False
    summary_tsv: Path | None = None
    metadata_json: Path | None = None
    parallel_contigs: bool = False


@dataclass
class WindowIndexSummary:
    cohort: str
    sample: str
    status: str
    windows_dir: str
    total_windows: int
    kept_windows: int
    dropped_n: int
    dropped_blacklist: int
    kept_bp: int
    dropped_bp: int
    kept_windows_with_blacklist: int
    kept_blacklist_bp: int
    low_coord_coverage_kept: int
    elapsed_sec: float
    message: str = ""


def run_window_index_generation(config: WindowIndexConfig) -> dict[str, object]:
    output_root = config.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    chroms = expand_chroms(config.chroms)
    samples = load_sample_ids(config.sample_list, config.samples)
    if not samples:
        raise ValueError("sample list is empty")
    policy = normalize_sequence_variant_policy(config.sequence_variant_policy)
    if config.sequence_source not in {"vcf_consensus", "assembly_direct"}:
        raise ValueError("--sequence-source must be vcf_consensus or assembly_direct")
    # Python API caller가 blacklist_bed=None으로 assembly 경로를 명시한 경우는
    # CLI의 --disable-reference-blacklist와 같은 의도다. 복제 config를 써서
    # 이후 metadata/reuse contract도 일관되게 "disabled"로 남긴다.
    if config.sequence_source == "assembly_direct" and config.blacklist_bed is None:
        config = replace(config, disable_reference_blacklist=True)
    if config.disable_reference_blacklist and config.sequence_source != "assembly_direct":
        raise ValueError("--disable-reference-blacklist is allowed only with --sequence-source assembly_direct")
    if config.sequence_source == "assembly_direct" and not config.disable_reference_blacklist:
        # assembly contig와 GRCh38의 multi-chain projection은 아직 제공하지 않는다.
        # 존재하지 않는 reference interval을 blacklist BED에 대조하지 않도록 source
        # 특이 정책 오류를 BED 파일 존재 검사보다 먼저 명확히 보고한다.
        raise ValueError("assembly_direct requires disabled GRCh38 blacklist until official multi-chain projection is implemented")
    if config.window_size <= 0:
        raise ValueError("window_size must be positive")
    offset_unit = effective_offset_unit(config)
    if offset_unit <= 0 or offset_unit > config.window_size:
        raise ValueError("offset_unit must be within 1..window_size")
    if config.max_blacklist_frac < 0:
        raise ValueError("max_blacklist_frac must be non-negative")

    config=replace(config,parallel_contigs=len(samples)==1)
    blacklist = {} if config.disable_reference_blacklist else load_blacklist(config.blacklist_bed)
    blacklist_md5 = "" if config.disable_reference_blacklist else (file_md5(config.blacklist_bed) if config.blacklist_bed else "")
    blacklist_bp = sum(end - start for values in blacklist.values() for start, end in values)
    started_at = time.time()
    summaries: list[WindowIndexSummary] = []
    with ThreadPoolExecutor(max_workers=max(1, config.jobs)) as executor:
        futures = {
            executor.submit(process_one_sample, config, sample, output_root, chroms, policy, blacklist, blacklist_md5, blacklist_bp): sample
            for sample in samples
        }
        for future in as_completed(futures):
            summaries.append(future.result())

    summaries.sort(key=lambda row: row.sample)
    summary_tsv = (config.summary_tsv or (output_root / f"{config.cohort}.window_index.summary.tsv")).resolve()
    metadata_json = (config.metadata_json or (output_root / f"{config.cohort}.window_index.metadata.json")).resolve()
    write_summary_tsv(summary_tsv, summaries)
    completed_at = time.time()
    failed = [row for row in summaries if row.status == "FAIL"]
    metadata = {
        "schema": WINDOW_INDEX_SCHEMA,
        "process_code_schema_compat": PROCESS_CODE_WINDOW_SCHEMA,
        "cohort": config.cohort,
        "sample_count": len(samples),
        "chroms": chroms,
        "output_root": str(output_root),
        "sequence_variant_policy": policy,
        "sequence_source": config.sequence_source,
        "reference_blacklist_enabled": not config.disable_reference_blacklist,
        "window_size": config.window_size,
        "offset_unit": offset_unit,
        "seed": config.seed,
        "window_start_policy": "deterministic_random_offset_sha256_seed_sample_haplotype_contig",
        "drop_n": bool(config.drop_n),
        "drop_reason_priority": DROP_REASON_PRIORITY,
        "n_fraction_denominator": N_FRACTION_DENOMINATOR,
        "blacklist_bed": str(config.blacklist_bed or ""),
        "blacklist_md5": blacklist_md5,
        "blacklist_regions": sum(len(values) for values in blacklist.values()),
        "blacklist_bp": blacklist_bp,
        "max_blacklist_frac": config.max_blacklist_frac,
        "coord_coverage_filter": None,
        "coord_coverage_policy": "record_only_not_filter",
        "summary_tsv": str(summary_tsv),
        "started_at_unix": started_at,
        "completed_at_unix": completed_at,
        "elapsed_sec": round(completed_at - started_at, 3),
        "ok_count": sum(1 for row in summaries if row.status == "OK"),
        "skip_count": sum(1 for row in summaries if row.status == "SKIP"),
        "fail_count": len(failed),
    }
    write_json(metadata_json, metadata)
    if failed:
        labels = ", ".join(f"{row.sample}:{row.message}" for row in failed[:5])
        raise RuntimeError(f"window index generation failed: {labels}")
    return metadata


def process_one_sample(
    config: WindowIndexConfig,
    sample: str,
    output_root: Path,
    chroms: list[str],
    policy: str,
    blacklist: dict[str, list[tuple[int, int]]],
    blacklist_md5: str,
    blacklist_bp: int,
) -> WindowIndexSummary:
    t0 = time.time()
    windows_dir = output_root / sample / "windows"
    tsv = windows_dir / "windows.tsv"
    blacklist_tsv = windows_dir / "windows_blacklist.tsv"
    meta_path = windows_dir / "windows.json"
    try:
        if tsv.exists() and blacklist_tsv.exists() and meta_path.exists() and not config.overwrite:
            try:
                previous = json.loads(meta_path.read_text(encoding="utf-8"))
            except Exception:
                previous = {}
            if reusable_metadata(previous, config, chroms, policy, blacklist_md5):
                # 좌표 열 검증 없이 구형 캐시를 재사용하면 RNA 단계에서 실패한다.
                # 행·kept·offset은 유지하고 동치인 personal 좌표 열만 보충한다.
                ensure_window_contract(tsv)
                ensure_window_n_annotations(tsv)
                canonical = n_value_metadata(previous)
                if canonical != previous:
                    previous = canonical
                    write_json(meta_path, previous)
                if previous.get("drop_reason_priority") != DROP_REASON_PRIORITY:
                    # 유지 창과 TSV 해시는 보존하고, 구형 캐시의 제외 통계만 갱신한다.
                    cached_rows = read_windows(tsv)
                    for hap, stats in previous.get("haplotypes", {}).items():
                        stats.update(count_index_exclusions(
                            [row for row in cached_rows if row["haplotype"] == hap],
                            max_blacklist_frac=config.max_blacklist_frac,
                            drop_n=config.drop_n,
                        ))
                    previous["drop_reason_priority"] = DROP_REASON_PRIORITY
                    previous["n_fraction_denominator"] = N_FRACTION_DENOMINATOR
                    write_json(meta_path, previous)
                hap_stats = previous.get("haplotypes", {}) if isinstance(previous, dict) else {}
                aggregate = aggregate_previous_stats(hap_stats)
                return WindowIndexSummary(
                    cohort=config.cohort,
                    sample=sample,
                    status="SKIP",
                    windows_dir=str(windows_dir),
                    elapsed_sec=round(time.time() - t0, 3),
                    message="existing compatible window index",
                    **aggregate,
                )
        seq_root = output_root / sample / "individual_fasta"
        if not seq_root.is_dir():
            raise FileNotFoundError(f"individual_fasta directory is missing: {seq_root}")
        windows_dir.mkdir(parents=True, exist_ok=True)
        # assembly_direct에는 GRCh38↔assembly multi-chain 투영이 없으므로 reference
        # blacklist와 chain-coordinate coverage를 쓰지 않는다. VCF consensus는 기존처럼
        # chain을 필수로 유지한다.
        chain_required = policy != "snp_only" and config.sequence_source != "assembly_direct"
        all_rows: list[dict[str, Any]] = []
        all_blacklist_rows: list[dict[str, Any]] = []
        hap_stats: dict[str, dict[str, int]] = {}
        tasks=[]
        for hap in ("HAP1", "HAP2"):
            hap_stats[hap]=empty_stats()
            for contig in chroms:
                tasks.append(dict(sample=sample,hap=hap,contig=contig,
                    fasta=find_fasta(seq_root/hap,contig),
                    chain=find_chain(seq_root/hap,contig,required=chain_required),
                    config=config,blacklist=blacklist))
        from .window_parallel import run as window_jobs
        for hap,rows,bl_rows,contig_stats in window_jobs(tasks,parallel=config.parallel_contigs):
            all_rows.extend(rows)
            all_blacklist_rows.extend(bl_rows)
            add_stats(hap_stats[hap],contig_stats)
        write_tsv(
            tsv,
            all_rows,
            [
                "contig",
                "haplotype",
                "window_index",
                "raw_start",
                "actual_length",
                # raw_*는 compact/trace legacy contract, personal_*는 RNA derived
                # coverage/final package contract다. VCF consensus HAP의 raw axis는
                # personal HAP axis이므로 같은 interval을 두 이름으로 명시한다.
                "personal_start",
                "personal_end",
                # 한 HAP window는 indel/chain gap 때문에 단일 REF interval로 항상
                # 표현되지 않는다. exact REF projection은 blacklist TSV에 segment로
                # 남기고, 여기서는 non-contiguous sentinel을 빈 값으로 둔다.
                "ref_start",
                "ref_end",
                *N_VALUE_FIELDS, "storage_policy", "storage_drop_reason",
                "coord_coverage",
                "bl_frac",
                "kept",
                "drop_reason",
            ],
        )
        write_tsv(
            blacklist_tsv,
            all_blacklist_rows,
            ["contig", "haplotype", "window_index", "personal_start", "personal_end", "ref_start", "ref_end"],
        )
        offset_unit = effective_offset_unit(config)
        meta = {
            "schema": WINDOW_INDEX_SCHEMA,
            "process_code_schema_compat": PROCESS_CODE_WINDOW_SCHEMA,
            "cohort": config.cohort,
            "sample": sample,
            "chroms": chroms,
            "sequence_variant_policy": policy,
            "sequence_source": config.sequence_source,
            "reference_blacklist_enabled": not config.disable_reference_blacklist,
            "window_bp": config.window_size,
            "overlap": "none",
            "coordinate_space": "individual_raw_0based_halfopen",
            "window_start_policy": "deterministic_random_offset_sha256_seed_sample_haplotype_contig",
            "offset_unit": offset_unit,
            "seed": config.seed,
            "drop_n": bool(config.drop_n),
            "drop_reason_priority": DROP_REASON_PRIORITY,
            "n_fraction_denominator": N_FRACTION_DENOMINATOR,
            "n_policy": "drop any non-ACGT" if config.drop_n else "preserve internal N; exclude windows without ACGT",
            "storage_policy": "annotated_acgtn_v1",
            "all_n_windows_excluded": True,
            "training_filter_applied_to_storage": False,
            "coord_coverage_filter": None,
            "coord_coverage_policy": "record only; not used for kept decision",
            "blacklist_bed": str(config.blacklist_bed or ""),
            "blacklist_md5": blacklist_md5,
            "blacklist_regions": sum(len(values) for values in blacklist.values()),
            "blacklist_bp": blacklist_bp,
            "max_blacklist_frac": config.max_blacklist_frac,
            "blacklist_policy": (
                "disabled_for_assembly_direct_no_multichain_projection"
                if config.disable_reference_blacklist
                else "reference-coordinate overlap; chain-lift individual windows for snv_indel; no masking; residual hits recorded"
            ),
            "residual_blacklist_file": "windows_blacklist.tsv",
            "windows_file": "windows.tsv",
            "haplotypes": hap_stats,
        }
        meta = n_value_metadata(meta)
        write_json(meta_path, meta)
        aggregate = aggregate_previous_stats(hap_stats)
        return WindowIndexSummary(
            cohort=config.cohort,
            sample=sample,
            status="OK",
            windows_dir=str(windows_dir),
            elapsed_sec=round(time.time() - t0, 3),
            message="created",
            **aggregate,
        )
    except Exception as exc:  # noqa: BLE001 - sample failures are reported in summary TSV.
        return WindowIndexSummary(
            cohort=config.cohort,
            sample=sample,
            status="FAIL",
            windows_dir=str(windows_dir),
            total_windows=0,
            kept_windows=0,
            dropped_n=0,
            dropped_blacklist=0,
            kept_bp=0,
            dropped_bp=0,
            kept_windows_with_blacklist=0,
            kept_blacklist_bp=0,
            low_coord_coverage_kept=0,
            elapsed_sec=round(time.time() - t0, 3),
            message=str(exc),
        )


def windows_for_fasta(
    *,
    sample: str,
    hap: str,
    contig: str,
    fasta: Path,
    chain: Path | None,
    config: WindowIndexConfig,
    blacklist: dict[str, list[tuple[int, int]]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, int]]:
    query_name, original = read_fasta_sequence(fasta, contig)
    trimmed = trim_terminal_n(original)
    sequence = trimmed.sequence.upper()
    offset = deterministic_offset(config.seed, sample, f"{hap}:{contig}", effective_offset_unit(config))
    alignment = load_chain_alignment(chain, query_name=query_name, target_name=contig) if chain else None
    mapped_intervals = chain_query_intervals(alignment, len(original))
    chain_blocks = chain_blacklist_blocks(alignment) if alignment else None
    rows: list[dict[str, Any]] = []
    bl_rows: list[dict[str, Any]] = []
    stats = empty_stats()
    bl_iv = blacklist.get(contig, [])
    window_index = 0
    for local_start in range(offset, len(sequence), config.window_size):
        actual = sequence[local_start: local_start + config.window_size]
        if not actual:
            continue
        raw_start = trimmed.source_start + local_start
        raw_end = raw_start + len(actual)
        annotations = sequence_annotations(actual)
        n_nonacgt = annotations["n_nonacgt"]
        mapped_bases, coverage = interval_mapping_coverage(mapped_intervals, raw_start, raw_end)
        bl_frac, segments = blacklist_hits(bl_iv, raw_start, raw_end, chain_blocks)
        ok_n = n_nonacgt < len(actual) and ((n_nonacgt == 0) or (not config.drop_n))
        ok_bl = bl_frac <= config.max_blacklist_frac
        kept = ok_n and ok_bl
        drop_reason = exclusion_reason(ok_n=ok_n, ok_blacklist=ok_bl)
        stats["total_windows"] += 1
        if kept:
            stats["kept_windows"] += 1
            stats["kept_bp"] += len(actual)
            if coverage < 0.95:
                stats["low_coord_coverage_kept"] += 1
            if segments:
                stats["kept_windows_with_blacklist"] += 1
                for personal_start, personal_end, ref_start, ref_end in segments:
                    stats["kept_blacklist_bp"] += personal_end - personal_start
                    bl_rows.append(
                        {
                            "contig": contig,
                            "haplotype": hap,
                            "window_index": window_index,
                            "personal_start": personal_start,
                            "personal_end": personal_end,
                            "ref_start": ref_start,
                            "ref_end": ref_end,
                        }
                    )
        else:
            stats["dropped_" + drop_reason] += 1
            stats["dropped_bp"] += len(actual)
        rows.append(
            {
                "contig": contig,
                "haplotype": hap,
                "window_index": window_index,
                "raw_start": raw_start,
                "actual_length": len(actual),
                "personal_start": raw_start,
                "personal_end": raw_end,
                "ref_start": "",
                "ref_end": "",
                **annotations,
                "storage_drop_reason": "all_n" if n_nonacgt == len(actual) else drop_reason,
                "coord_coverage": f"{coverage:.6f}",
                "bl_frac": f"{bl_frac:.6f}",
                "kept": int(kept),
                "drop_reason": drop_reason,
            }
        )
        _ = mapped_bases  # retained for readability; coverage is the stored field.
        window_index += 1
    return rows, bl_rows, stats


def read_fasta_sequence(path: Path, contig: str) -> tuple[str, str]:
    if not path.exists() or not path.stat().st_size:
        raise FileNotFoundError(f"FASTA is missing: {path}")
    fasta = pysam.FastaFile(str(path))
    try:
        candidates = [contig, contig.removeprefix("chr")]
        for name in candidates:
            if name in fasta.references:
                return name, fasta.fetch(name).upper()
        if len(fasta.references) == 1:
            name = fasta.references[0]
            return name, fasta.fetch(name).upper()
        raise KeyError(f"contig {contig} not found in {path}; references={fasta.references[:5]}")
    finally:
        fasta.close()


def load_blacklist(path: Path | None) -> dict[str, list[tuple[int, int]]]:
    if path is None:
        return {}
    if not path.exists() or not path.stat().st_size:
        raise FileNotFoundError(f"blacklist BED is missing: {path}")
    out: dict[str, list[tuple[int, int]]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            out.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    for values in out.values():
        values.sort()
    return out


def chain_query_intervals(alignment: ChainAlignment | None, raw_length: int) -> list[tuple[int, int]]:
    if alignment is None:
        return [(0, raw_length)]
    intervals = [(int(q), int(q) + int(length)) for q, _t, length in alignment.blocks]
    return merge_intervals(intervals)


def chain_blacklist_blocks(alignment: ChainAlignment | None) -> list[tuple[int, int, int]] | None:
    if alignment is None:
        return None
    return [(int(q), int(q) + int(length), int(t)) for q, t, length in alignment.blocks]


def blacklist_hits(
    intervals: list[tuple[int, int]],
    start: int,
    end: int,
    blocks: list[tuple[int, int, int]] | None,
) -> tuple[float, list[tuple[int, int, int, int]]]:
    if end <= start or not intervals:
        return 0.0, []
    hit = 0
    segments: list[tuple[int, int, int, int]] = []
    if blocks is None:
        for ref_start, ref_end in overlap_segments(intervals, start, end):
            hit += ref_end - ref_start
            segments.append((ref_start, ref_end, ref_start, ref_end))
        return hit / (end - start), segments
    for query_start, query_end, target_start in blocks:
        if query_end <= start:
            continue
        if query_start >= end:
            break
        personal_start = max(start, query_start)
        personal_end = min(end, query_end)
        ref_start = target_start + (personal_start - query_start)
        ref_end = target_start + (personal_end - query_start)
        for overlap_start, overlap_end in overlap_segments(intervals, ref_start, ref_end):
            hit += overlap_end - overlap_start
            segments.append(
                (
                    personal_start + (overlap_start - ref_start),
                    personal_start + (overlap_end - ref_start),
                    overlap_start,
                    overlap_end,
                )
            )
    segments.sort()
    return hit / (end - start), segments


def overlap_segments(intervals: list[tuple[int, int]], start: int, end: int) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for left, right in intervals:
        if right <= start:
            continue
        if left >= end:
            break
        overlap_start = max(start, left)
        overlap_end = min(end, right)
        if overlap_start < overlap_end:
            out.append((overlap_start, overlap_end))
    return out


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    cleaned = sorted((int(start), int(end)) for start, end in intervals if int(end) > int(start))
    if not cleaned:
        return []
    merged = [cleaned[0]]
    for start, end in cleaned[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def find_fasta(directory: Path, contig: str) -> Path:
    candidates = [directory / f"{contig}{suffix}" for suffix in (".fa", ".fasta", ".fa.gz", ".fasta.gz")]
    candidates.extend(directory / contig / f"{contig}{suffix}" for suffix in (".fa", ".fasta", ".fa.gz", ".fasta.gz"))
    for path in candidates:
        if path.exists():
            return path
    shown = ", ".join(str(path) for path in candidates[:4])
    raise FileNotFoundError(f"missing FASTA for {contig} under {directory}; first candidates: {shown}")


def find_chain(directory: Path, contig: str, *, required: bool) -> Path | None:
    candidates = [directory / f"{contig}.chain", directory / contig / f"{contig}.chain"]
    for path in candidates:
        if path.exists() and path.stat().st_size:
            return path
    if required:
        shown = ", ".join(str(path) for path in candidates)
        raise FileNotFoundError(f"missing chain for {contig} under {directory}; candidates: {shown}")
    return None


def reusable_metadata(previous: dict[str, Any], config: WindowIndexConfig, chroms: list[str], policy: str, blacklist_md5: str) -> bool:
    return (
        previous.get("schema") == WINDOW_INDEX_SCHEMA
        and previous.get("storage_policy") == "annotated_acgtn_v1"
        and previous.get("chroms") == chroms
        and previous.get("sequence_variant_policy") == policy
        and previous.get("sequence_source", "vcf_consensus") == config.sequence_source
        and bool(previous.get("reference_blacklist_enabled", True)) == (not config.disable_reference_blacklist)
        and int(previous.get("window_bp", -1)) == int(config.window_size)
        and int(previous.get("offset_unit", -1)) == int(effective_offset_unit(config))
        and int(previous.get("seed", -1)) == int(config.seed)
        and bool(previous.get("drop_n")) == bool(config.drop_n)
        and float(previous.get("max_blacklist_frac", -1.0)) == float(config.max_blacklist_frac)
        and str(previous.get("blacklist_md5", "")) == str(blacklist_md5)
    )


def empty_stats() -> dict[str, int]:
    return {
        "total_windows": 0,
        "kept_windows": 0,
        "dropped_n": 0,
        "dropped_blacklist": 0,
        "kept_bp": 0,
        "dropped_bp": 0,
        "kept_windows_with_blacklist": 0,
        "kept_blacklist_bp": 0,
        "low_coord_coverage_kept": 0,
    }


def add_stats(target: dict[str, int], source: dict[str, int]) -> None:
    for key in target:
        target[key] += int(source.get(key, 0))


def aggregate_previous_stats(hap_stats: dict[str, Any]) -> dict[str, int]:
    total = empty_stats()
    for value in hap_stats.values():
        if isinstance(value, dict):
            add_stats(total, {key: int(value.get(key, 0)) for key in total})
    return total


def effective_offset_unit(config: WindowIndexConfig) -> int:
    return int(config.window_size if config.offset_unit is None else config.offset_unit)


def file_md5(path: Path | None) -> str:
    if path is None:
        return ""
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_summary_tsv(path: Path, rows: list[WindowIndexSummary]) -> None:
    write_tsv(path, [asdict(row) for row in rows], list(WindowIndexSummary.__dataclass_fields__.keys()))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
