#!/usr/bin/env python3
"""WGS Stage 23/RNA Stage 25: cleanup 전에 sample package를 완전 검증한다."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import zarr

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE_ROOT))
sys.path.insert(0, str(PIPELINE_ROOT / "src"))
from integrate_wgs_rna.window_contract import read_windows
from integrate_wgs_rna.validation_parallel import hash_workers, sequence_checks, depth_checks, uint32_stats as zarr_uint32_stats
from tools.compact_manifest import completed_contigs
from tools.personal_assets import MANIFEST as PERSONAL_ASSET_MANIFEST
from tools.personal_assets import SCHEMA as PERSONAL_ASSET_SCHEMA
from tools.personal_assets import build as build_personal_assets
from tools.personal_assets import verify as verify_personal_assets


def sha256(path: Path, chunk: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    with path.open() as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"JSON object가 아님: {path}")
    return value


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(f".{path.name}.part.{os.getpid()}")
    part.write_text(text)
    os.replace(part, path)


class Audit:
    def __init__(self) -> None:
        self.checks: list[dict] = []

    def check(self, name: str, ok: bool, detail: str = "") -> None:
        self.checks.append({"name": name, "status": "pass" if ok else "fail", "detail": detail})

    @property
    def passed(self) -> bool:
        return all(item["status"] == "pass" for item in self.checks)


def validate_preflight(path: Path, audit: Audit, name: str, sample: str, schema: str) -> None:
    try:
        data = read_json(path)
        ok = data.get("schema") == schema and data.get("sample") == sample and data.get("status") == "pass"
        audit.check(name, ok, f"schema={data.get('schema')} sample={data.get('sample')} status={data.get('status')}")
    except Exception as exc:
        audit.check(name, False, str(exc))


def validate_personal_assets(root: Path, audit: Audit, contigs: set[str]) -> str | None:
    """cleanup까지 보존할 HAP FASTA/FAI/chain/QC를 final PASS 전에 검사한다."""
    try:
        asset_root = root / "individual_fasta"
        manifest_path = asset_root / PERSONAL_ASSET_MANIFEST
        expected_files = 1 + 2 * len(contigs) * 3
        with hash_workers() as workers:
            if manifest_path.is_file():
                manifest = read_json(manifest_path)
                verified = verify_personal_assets(asset_root,workers=workers)
                result = {
                    **manifest,
                    "status": verified.get("status"),
                    "errors": verified.get("problems", []),
                }
            else:
                # 구형 output의 1회 migration: 구조·chain/reference를 모두 검사한 뒤에만
                # 현재 자산을 기준점으로 고정한다. 신규 output은 Stage 09에서 이미 생성된다.
                result = build_personal_assets(asset_root, sorted(contigs), None,workers=workers)
                if result.get("status") == "pass":
                    atomic_text(
                        manifest_path,
                        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                    )
        manifest_contigs = {str(value) for value in result.get("contigs", [])}
        ok = (
            result.get("schema") == PERSONAL_ASSET_SCHEMA
            and result.get("status") == "pass"
            and manifest_contigs == contigs
            and len(result.get("files", {})) == expected_files
        )
        audit.check(
            "personal_hap_assets", ok,
            f"schema={result.get('schema')} status={result.get('status')} "
            f"contigs={len(manifest_contigs)}/{len(contigs)} "
            f"files={len(result.get('files', {}))}/{expected_files} "
            f"errors={result.get('errors', [])[:3]}",
        )
        return sha256(manifest_path) if ok else None
    except Exception as exc:
        audit.check("personal_hap_assets", False, str(exc))
        return None


def validate_provenance(path: Path, audit: Audit, sample: str, mode: str, root: Path) -> None:
    """최종 validator 이전 단계가 모두 success/conditional-skip로 끝났는지 확인한다."""
    original_layout = (PIPELINE_ROOT / "09_personal_fasta.sh").is_file()
    wgs_stop, rna_start, rna_stop = (13, 13, 25) if original_layout else (11, 11, 23)
    expected = {"WGS": {f"{stage:02d}" for stage in range(1, wgs_stop)}}
    if mode == "rna":
        expected["RNA"] = {f"{stage:02d}" for stage in range(rna_start, rna_stop)}
    completed = {branch: set() for branch in expected}
    schemas_ok = True
    log_errors: list[str] = []
    records = 0
    try:
        with path.open() as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                records += 1
                data = json.loads(line)
                if data.get("schema") != "mfm_execution_provenance_v1" or data.get("sample_id") != sample:
                    schemas_ok = False
                    continue
                branch, stage, event = str(data.get("branch", "")), str(data.get("stage", "")), data.get("event")
                if branch not in expected or stage not in expected[branch]:
                    continue
                if event == "finish" and data.get("status") == "success" and int(data.get("exit_code", -1)) == 0:
                    completed[branch].add(stage)
                    log_path = Path(str(data.get("log_path", "")))
                    if not log_path.is_file() or root.resolve() / "logs" / "stage_logs" not in log_path.resolve().parents:
                        log_errors.append(f"{branch}:{stage}:{log_path}")
                elif event == "skip" and data.get("status") == "skipped":
                    completed[branch].add(stage)
    except Exception as exc:
        audit.check("execution_provenance", False, str(exc))
        return
    missing = {branch: sorted(stages - completed[branch]) for branch, stages in expected.items() if stages - completed[branch]}
    audit.check(
        "execution_provenance",
        records > 0 and schemas_ok and not missing,
        f"records={records} missing={missing} schema_sample_ok={schemas_ok}",
    )
    audit.check("stage_logs", not log_errors, f"missing_or_external={log_errors[:5]}")


def load_windows(path: Path, expected_contigs: set[str]) -> tuple[list[dict], dict[str, list[str]]]:
    rows = read_windows(path)
    keys: dict[str, list[str]] = {"HAP1": [], "HAP2": []}
    kept_rows: list[dict] = []
    observed: dict[str, set[str]] = {"HAP1": set(), "HAP2": set()}
    for row in rows:
        if str(row["kept"]).strip().lower() not in {"1", "true"}:
            continue
        kept_rows.append(row)
        hap = row["haplotype"]
        if hap not in keys or row["contig"] not in expected_contigs:
            raise ValueError(f"kept window에 예상하지 않은 HAP/contig: {row}")
        start, end = int(row["personal_start"]), int(row["personal_end"])
        if start < 0 or end <= start:
            raise ValueError(f"잘못된 window 좌표: {row}")
        key = f'{row["contig"]}:{hap}:{int(row["window_index"])}:{start}-{end}'
        keys[hap].append(key)
        observed[hap].add(row["contig"])
    for hap in keys:
        if observed[hap] != expected_contigs:
            raise ValueError(
                f"{hap} window contig 불일치 expected={sorted(expected_contigs)} observed={sorted(observed[hap])}"
            )
        if len(keys[hap]) != len(set(keys[hap])):
            raise ValueError(f"{hap} window key 중복")
    if not kept_rows:
        raise ValueError("kept window가 없음")
    return rows, keys


def validate_compact(a: argparse.Namespace, audit: Audit, contigs: set[str]) -> Path:
    layout = a.root / "compact" / a.layout
    manifest_path = layout / "manifest.json"
    audit.check("compact_layout", layout.is_dir(), str(layout))
    if not manifest_path.is_file():
        audit.check("compact_manifest", False, str(manifest_path))
        return layout
    try:
        manifest = read_json(manifest_path)
        payload=getattr(a,'wgs_payload','full')
        if manifest.get('wgs_payload','full')!=payload:raise ValueError('요청한 WGS payload와 산출물 불일치')
        if payload=='pileup_feature':
            from integrate_wgs_rna.pileup_features.schema import SCHEMA, LEGACY_SCHEMA, LEGACY_CHANNELS, resolve_channels
            if manifest.get('schema')==SCHEMA:
                resolve_channels(manifest.get('channels'),manifest.get('feature_set'))
                a.pileup_feature_set=manifest['feature_set']
            elif manifest.get('schema')==LEGACY_SCHEMA and manifest.get('channels')==list(LEGACY_CHANNELS):
                a.pileup_feature_set='legacy_reference12_v1'
            else:raise ValueError('알 수 없는 pileup 채널 계약')
            metadata=read_json(a.root/'metadata/sample.json')
            if manifest.get('schema')==SCHEMA and (metadata.get('pileup_feature_set')!=a.pileup_feature_set or metadata.get('pileup_feature_schema')!=SCHEMA):
                raise ValueError('pileup metadata/실제 채널 계약 불일치')
        passed = completed_contigs(manifest, layout, a.sample)
        ok = (
            manifest.get("status") == "complete"
            and manifest.get("validation_level") == "full"
            and passed == contigs
        )
        audit.check(
            "compact_manifest", ok,
            f"status={manifest.get('status')} validation={manifest.get('validation_level')} contigs={len(passed)}",
        )
        global_rows = read_windows(a.root / "windows" / "windows.tsv")
        ordered_contigs=sorted(contigs)
        sequence_results=(sequence_checks([layout/c for c in ordered_contigs],source=payload=='sequence_only')
                          if payload in {'sequence_only','pileup_feature'} else [])
        for index,contig in enumerate(ordered_contigs):
            directory = layout / contig
            required = [
                "contig_metadata.json", "windows.tsv", "hap_sequence.zarr", "depth.zarr",
                "readset.HAP1.zarr", "readset.HAP2.zarr",
            ]
            if payload=='sequence_only':
                required=['contig_metadata.json','windows.tsv','hap_sequence.zarr','filter_qc.json']
                if sequence_results[index][1] is not None:raise ValueError(sequence_results[index][1])
            if payload=='pileup_feature':
                required=['contig_metadata.json','windows.tsv','hap_sequence.zarr','pileup.zarr','filter_qc.json']
                if sequence_results[index][1] is not None:raise ValueError(sequence_results[index][1])
            missing = [name for name in required if not (directory / name).exists()]
            audit.check(f"compact:{contig}", not missing, ",".join(missing))
            # compact index는 kept window를 재번호화한다. 위치·길이·행 순서만
            # 비교해야 전역 audit index와의 번호 차이를 corruption으로 오판하지 않는다.
            compact_rows = read_windows(directory / "windows.tsv", compact=True)
            for hap in ("HAP1", "HAP2"):
                def spans(rows):
                    return [(int(r["personal_start"]), int(r["personal_end"])) for r in rows
                            if r["contig"] == contig and r["haplotype"] == hap
                            and r["kept"].strip().lower() in {"1", "true"}]
                audit.check(f"compact:{contig}:{hap}:window_contract",
                            spans(global_rows) == spans(compact_rows), "개인 좌표·길이·순서 대응")
    except Exception as exc:  # 실패 보고서를 남기기 위해 예외를 check로 변환한다.
        audit.check("compact_manifest", False, str(exc))
    return layout


def validate_trace(a: argparse.Namespace, audit: Audit, contigs: set[str], layout: Path) -> None:
    if a.trace_policy == "off" and a.allow_without_trace:
        audit.check("source_trace", True, "explicit override: trace policy off")
        return
    trace_mode=a.wgs_payload if getattr(a,'wgs_payload','full') in {'sequence_only','pileup_feature'} else a.hap_read_mode
    path = a.root / "qc" / "stage12_trace" / trace_mode / "trace_manifest.json"
    if not path.is_file():
        audit.check("source_trace", False, str(path))
        return
    try:
        data = read_json(path)
        got = {key for key, value in data.get("contigs", {}).items() if value == "pass"}
        root_matches = Path(data.get("compact_root", "/")).resolve() == layout.resolve()
        audit.check(
            "source_trace",
            data.get("status") == "pass" and got == contigs and root_matches,
            f"status={data.get('status')} contigs={len(got)} root_matches={root_matches}",
        )
        newer = next((p for p in layout.rglob("*") if p.is_file() and p.stat().st_mtime_ns > path.stat().st_mtime_ns), None)
        audit.check("source_trace_fresh", newer is None, str(newer or ""))
    except Exception as exc:
        audit.check("source_trace", False, str(exc))


def validate_rna(a: argparse.Namespace, audit: Audit, window_keys: dict[str, list[str]]) -> None:
    """RNA 전장 evidence와 WGS-coordinate derived coverage를 각각 검증한다."""
    rna_root = a.root / "rna_personal" / a.rna_label
    if getattr(a,'rna_payload','full')=='coverage_only':
        try:
            from rna_signal import validate_signal
            validate_signal(a.root,a.sample,a.rna_label)
            audit.check('rna_signal_evidence_contract',True,'RNA signal 전용: ASE/CNV 비적용')
        except Exception as exc:
            audit.check('rna_signal_evidence_contract',False,str(exc))
    else:
        required = [
            rna_root / "hap_variants.tsv",
            rna_root / "hap_variants_qc.json",
            rna_root / "allelic_counts.tsv",
            rna_root / "allelic_counts.meta.json",
            rna_root / "allelic_counts.split_half.tsv",
            rna_root / "allelic_counts.both_haps.tsv",
            rna_root / "allelic_counts.mmap_adjusted.tsv",
            rna_root / "ase_qc.json",
            rna_root / "ase_training_eligibility.tsv",
            rna_root / "rna_realign_summary.json",
        ]
        for path in required:
            audit.check(f"rna_output:{path.name}", path.is_file() and path.stat().st_size > 0, str(path))

        # Stage 20 → 21 → 23의 근거 파일이 실제로 같은 내용인지 hash/schema로 연결한다.
        hv_path = rna_root / "hap_variants.tsv"
        hv_qc_path = rna_root / "hap_variants_qc.json"
        count_meta_path = rna_root / "allelic_counts.meta.json"
        ase_qc_path = rna_root / "ase_qc.json"
        if all(path.is_file() for path in (hv_path, hv_qc_path, count_meta_path, ase_qc_path)):
            try:
                hv_sha = sha256(hv_path)
                hv_qc = read_json(hv_qc_path)
                count_meta = read_json(count_meta_path)
                ase_qc = read_json(ase_qc_path)
                count_paths = {
                    "counts": rna_root / "allelic_counts.tsv",
                    "split_half": rna_root / "allelic_counts.split_half.tsv",
                    "both_haps": rna_root / "allelic_counts.both_haps.tsv",
                    "mmap_adjusted": rna_root / "allelic_counts.mmap_adjusted.tsv",
                }
                observed_count_hashes = {name: sha256(path) for name, path in count_paths.items()}
                expected_count_hashes = count_meta.get("output_sha256", {})
                audit.check(
                    "rna_hap_variant_qc_contract",
                    hv_qc.get("schema") == "mfm_hap_variants_qc_v1"
                    and hv_qc.get("output_sha256") == hv_sha,
                    f"schema={hv_qc.get('schema')} hash_match={hv_qc.get('output_sha256') == hv_sha}",
                )
                audit.check(
                    "rna_allelic_v5_contract",
                    count_meta.get("schema") == "mfm_allelic_counts_meta_v5"
                    and count_meta.get("script_revision") == "2026-10-02_validated_hap_variant_evidence"
                    and count_meta.get("hap_variants_sha256") == hv_sha,
                    f"schema={count_meta.get('schema')} revision={count_meta.get('script_revision')} "
                    f"hash_match={count_meta.get('hap_variants_sha256') == hv_sha}",
                )
                sample_qc = ase_qc.get("sample_qc", {})
                ase_inputs = ase_qc.get("inputs", {})
                meta_sha = sha256(count_meta_path)
                input_hash_ok = (
                    expected_count_hashes == observed_count_hashes
                    and ase_inputs.get("counts_sha256") == observed_count_hashes["counts"]
                    and ase_inputs.get("split_half_sha256") == observed_count_hashes["split_half"]
                    and ase_inputs.get("both_haps_sha256") == observed_count_hashes["both_haps"]
                    and ase_inputs.get("mmap_adjusted_sha256") == observed_count_hashes["mmap_adjusted"]
                    and ase_inputs.get("meta_sha256") == meta_sha
                )
                cnv_path_raw = ase_inputs.get("cnv_rd")
                if cnv_path_raw:
                    cnv_path = Path(cnv_path_raw)
                    cnv_meta_path = Path(str(cnv_path) + ".meta.json")
                    cnv_meta = read_json(cnv_meta_path)
                    input_hash_ok = (
                        input_hash_ok
                        and ase_inputs.get("cnv_rd_sha256") == sha256(cnv_path)
                        and ase_inputs.get("cnv_meta_sha256") == sha256(cnv_meta_path)
                        and cnv_meta.get("output_sha256") == sha256(cnv_path)
                    )
                audit.check(
                    "rna_ase_validated_evidence_contract",
                    ase_qc.get("schema") == "mfm_ase_qc_v3"
                    and ase_qc.get("script_revision")
                        == "2026-10-02_validated_hap_variant_evidence_cnv_hap_consistency"
                    and ase_inputs.get("hap_variants_sha256") == hv_sha
                    and input_hash_ok
                    and sample_qc.get("provenance_ok") is True
                    and sample_qc.get("duplicates_excluded") is True,
                    f"schema={ase_qc.get('schema')} revision={ase_qc.get('script_revision')} "
                    f"hash_match={ase_inputs.get('hap_variants_sha256') == hv_sha} "
                    f"input_hashes={input_hash_ok} "
                    f"provenance={sample_qc.get('provenance_ok')} duplicates_excluded={sample_qc.get('duplicates_excluded')}",
                )
            except Exception as exc:
                audit.check("rna_validated_evidence_contract", False, str(exc))
    windows_path = a.root / "windows" / "windows.tsv"
    windows_sha = sha256(windows_path)
    window_rows = read_windows(windows_path)
    for hap in ("HAP1", "HAP2"):
        directory = rna_root / hap
        for name in ("annotation.gtf", "exons.bygene.tsv", "gtf_lift_qc.json", "rna_realign.json"):
            path = directory / name
            audit.check(f"rna_{hap}:{name}", path.is_file() and path.stat().st_size > 0, str(path))

        # 기존 v3 전장 track은 RNA QC/evidence용으로 유지한다.
        genome_store = directory / "coverage.zarr"
        genome_attrs = {}
        if not genome_store.is_dir():
            audit.check(f"rna_{hap}:genome_coverage", False, str(genome_store))
        else:
            try:
                genome = zarr.open_group(str(genome_store), mode="r")
                genome_attrs = dict(genome.attrs)
                requested_contigs = [contig for contig in a.chroms.replace(",", " ").split() if contig]
                present = all(contig in genome for contig in requested_contigs)
                identity_ok = (
                    str(genome.attrs.get("sample", a.sample)) == a.sample
                    and str(genome.attrs.get("hap", hap)) == hap
                )
                expected_digests = genome.attrs.get("depth_sha256_by_contig", {})
                observed_digests = {}
                total_depth = nonzero_positions = 0
                if present:
                    for contig, (total,nonzero,digest) in zip(requested_contigs,depth_checks(genome_store,requested_contigs)):
                        total_depth += total
                        nonzero_positions += nonzero
                        observed_digests[contig] = digest
                content_ok = (
                    expected_digests == observed_digests
                    and int(genome.attrs.get("total_depth", -1)) == total_depth
                    and int(genome.attrs.get("nonzero_positions", -1)) == nonzero_positions
                )
                audit.check(f"rna_{hap}:genome_coverage", present and identity_ok and content_ok,
                            f"contigs={len(list(genome.array_keys()))} identity={identity_ok} content={content_ok}")
            except Exception as exc:
                audit.check(f"rna_{hap}:genome_coverage", False, str(exc))

        store = directory / "coverage_windows.zarr"
        if not store.is_dir():
            audit.check(f"rna_{hap}:window_contract", False, str(store))
            continue
        try:
            group = zarr.open_group(str(store), mode="r")
            rows = [row for row in window_rows if row["haplotype"] == hap
                    and row["kept"].strip().lower() in {"1", "true"}]
            starts = [int(row["personal_start"]) for row in rows]
            ends = [int(row["personal_end"]) for row in rows]
            indices = [int(row["window_index"]) for row in rows]
            contig_names = [row["contig"] for row in rows]
            expected_offsets = [0]
            for start, end in zip(starts, ends):
                expected_offsets.append(expected_offsets[-1] + end - start)
            observed_keys = [str(value) for value in group["window_key"][:].tolist()]
            offsets = group["offsets"][:].tolist()
            window_total, window_nonzero, window_digest = zarr_uint32_stats(group["depth"])
            policy_keys = (
                "coverage_policy_revision", "duplicate_policy", "samtools_depth_args",
                "source_bam_signature", "source_fai_signature",
            )
            ok = (
                group.attrs.get("schema") == "mfm_rna_window_coverage_v1"
                and observed_keys == window_keys[hap]
                and offsets == expected_offsets
                and int(group["depth"].shape[0]) == expected_offsets[-1]
                and group["personal_start"][:].tolist() == starts
                and group["personal_end"][:].tolist() == ends
                and group["window_index"][:].tolist() == indices
                and [str(value) for value in group["contig"][:].tolist()] == contig_names
                and [str(value) for value in group["haplotype"][:].tolist()] == [hap] * len(rows)
                and group.attrs.get("windows_tsv_sha256") == windows_sha
                and str(group.attrs.get("sample")) == a.sample
                and str(group.attrs.get("hap")) == hap
                and all(group.attrs.get(key) == genome_attrs.get(key) for key in policy_keys)
                and int(group.attrs.get("total_depth", -1)) == window_total
                and int(group.attrs.get("nonzero_positions", -1)) == window_nonzero
                and group.attrs.get("depth_sha256") == window_digest
            )
            audit.check(f"rna_{hap}:window_contract", ok,
                        f"schema={group.attrs.get('schema')} windows={len(observed_keys)} "
                        f"depth={group['depth'].shape[0]}")
        except Exception as exc:
            audit.check(f"rna_{hap}:window_contract", False, str(exc))

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--mode", required=True, choices=("wgs", "rna"))
    parser.add_argument("--layout", required=True)
    parser.add_argument("--chroms", required=True)
    parser.add_argument("--hap-read-mode", required=True)
    parser.add_argument("--rna-payload", choices=("full","coverage_only"), default="full")
    parser.add_argument("--wgs-payload", choices=("full","sequence_only","pileup_feature"), default="full")
    parser.add_argument("--rna-label", default="")
    parser.add_argument("--trace-policy", default="required")
    parser.add_argument("--allow-without-trace", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--completeness-tsv", required=True, type=Path)
    args = parser.parse_args()

    audit = Audit()
    contigs = {token for token in args.chroms.replace(",", " ").split() if token}
    audit.check("sample_root", args.root.is_dir(), str(args.root))
    audit.check("sample_directory_name", args.root.name == args.sample, args.root.name)

    metadata = {}
    metadata_path = args.root / "metadata" / "sample.json"
    if metadata_path.is_file():
        try:
            metadata = read_json(metadata_path)
            observed = str(metadata.get("sample_id") or metadata.get("sample") or "")
            audit.check("metadata_sample_id", observed == args.sample, observed)
        except Exception as exc:
            audit.check("metadata_sample_id", False, str(exc))
    else:
        audit.check("metadata_sample_id", False, str(metadata_path))

    if metadata.get('rna_processing_status')=='deferred':
        audit.check('rna_deferred_contract',args.mode=='wgs' and args.wgs_payload in {'sequence_only','pileup_feature'}
                    and metadata.get('data_group')=='wgs_rna_pair' and bool(metadata.get('rna_path')),
                    'DNA 서열만 완료함. RNA coverage/ASE/CNV는 미처리이며 pair 입력 정보를 보존함.')
    windows_path = args.root / "windows" / "windows.tsv"
    window_keys: dict[str, list[str]] = {"HAP1": [], "HAP2": []}
    try:
        rows, window_keys = load_windows(windows_path, contigs)
        audit.check("wgs_windows", True, f"rows={len(rows)} sha256={sha256(windows_path)}")
    except Exception as exc:
        audit.check("wgs_windows", False, str(exc))

    window_metadata = args.root / "windows" / "windows.json"
    if window_metadata.exists() and read_json(window_metadata).get("storage_policy") == "annotated_acgtn_v1":
        try:
            from integrate_wgs_rna.window_policy import validate_window_annotations
            annotation_check = validate_window_annotations(args.root)
            audit.check("window_annotations_and_mapping", True, str(annotation_check))
        except Exception as exc:
            audit.check("window_annotations_and_mapping", False, str(exc))

    layout = validate_compact(args, audit, contigs)
    validate_trace(args, audit, contigs, layout)
    personal_assets_manifest_sha256 = validate_personal_assets(args.root, audit, contigs)
    preflight = args.root / "qc" / "wgs_preflight.json"
    if args.wgs_payload in {'sequence_only','pileup_feature'}:
        validate_preflight(args.root/'qc/sequence_preflight.json',audit,'sequence_preflight_qc',args.sample,'mfm_sequence_preflight_v1')
    if args.wgs_payload!='sequence_only':
        validate_preflight(preflight, audit, "wgs_preflight_qc", args.sample, "mfm_preflight_qc_v2")
    provenance = args.root / "metadata" / "execution_provenance.jsonl"
    validate_provenance(provenance, audit, args.sample, args.mode, args.root)

    if args.mode == 'rna':
        audit.check('rna_payload_contract',metadata.get('rna_payload','full')==args.rna_payload,str(metadata.get('rna_payload','full')))
        rna_preflight = args.root / "qc" / "rna_preflight.json"
        validate_preflight(
            rna_preflight, audit, "rna_preflight_qc", args.sample, "mfm_rna_preflight_contract_v2"
        )
        if not args.rna_label:
            audit.check("rna_label", False, "RNA mode에는 label 필요")
        else:
            validate_rna(args, audit, window_keys)

    status = "pass" if audit.passed else "fail"
    result = {
        "schema": "mfm_final_package_validation_v1", "sample": args.sample,
        "contract_revision": "2026-10-03_integrity_v2",
        "mode": args.mode, "status": status, "rna_payload": args.rna_payload,
        "validated_at": datetime.now(timezone.utc).isoformat(),
        "sample_root": str(args.root.resolve()), "chroms": sorted(contigs),
        "layout": args.layout, "hap_read_mode": args.hap_read_mode, "wgs_payload": args.wgs_payload,
        "personal_assets_manifest_sha256": personal_assets_manifest_sha256,
        "checks": audit.checks,
        "failed_checks": [item["name"] for item in audit.checks if item["status"] == "fail"],
    }
    if args.wgs_payload=='pileup_feature':result['pileup_feature_set']=getattr(args,'pileup_feature_set',None)
    if metadata.get('rna_processing_status')=='deferred':result['completion_scope']=metadata['completion_scope']
    if args.mode=='rna' and args.rna_payload=='coverage_only':
        signal=args.root/'rna_personal'/args.rna_label/'rna_signal_manifest.json'
        result['rna_signal_manifest_sha256']=sha256(signal) if signal.is_file() else None
        result['completion_scope']={'sequence_only':'dna_sequence_rna_coverage_only','pileup_feature':'dna_pileup_feature_rna_coverage_only','full':'wgs_full_rna_coverage_only'}[args.wgs_payload]
    atomic_text(args.output, json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    lines = ["sample\tmode\tcheck\tstatus\tdetail"]
    lines.extend(
        f'{args.sample}\t{args.mode}\t{item["name"]}\t{item["status"]}\t{item["detail"].replace(chr(9), " ")}'
        for item in audit.checks
    )
    atomic_text(args.completeness_tsv, "\n".join(lines) + "\n")
    print(json.dumps({"sample": args.sample, "mode": args.mode, "status": status,
                      "checks": len(audit.checks), "failed": result["failed_checks"]}, ensure_ascii=False))
    return 0 if audit.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
