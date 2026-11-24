"""최종 입력의 N/blacklist 주석과 개인→reference 매핑을 보존한다."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import shutil
from pathlib import Path

from .window_contract import (
    N_ANNOTATION_SCHEMA, LEGACY_N_FLAG_FIELDS, N_VALUE_FIELDS, N_VALUE_DEFINITIONS,
    ensure_window_n_annotations, n_value_metadata, validate_n_values,
)

POLICY = "annotated_acgtn_v1"


def _write_json(path: Path, value: dict) -> None:
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
    os.replace(part, path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def persist_window_annotations(sample_root: Path, compact_root: Path, index_dir: Path | None,
                               records: list[dict], reference_contract: dict | None,
                               sequence_source: str = "vcf_consensus") -> dict:
    """서열 배열을 다시 읽지 않고 후보 원장·주석·chain을 최종 패키지에 저장한다."""
    if index_dir is None:
        return {}
    source = index_dir / "windows.tsv"
    meta_path = index_dir / "windows.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    if "storage_policy" not in meta:
        return {"storage_policy": "legacy_external_index", "annotations_available": False}
    if meta.get("storage_policy") != POLICY:
        raise ValueError("새 compact 기본 정책에는 annotated_acgtn_v1 window index가 필요합니다")
    ensure_window_n_annotations(source)
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    metadata = sample_root / "metadata"
    metadata.mkdir(parents=True, exist_ok=True)
    for src, name in ((source, "wgs_windows.tsv"), (meta_path, "wgs_windows.json"),
                      (index_dir / "windows_blacklist.tsv", "wgs_windows_blacklist.tsv")):
        if src.exists():
            destination = metadata / name
            if src.resolve() != destination.resolve():
                shutil.copyfile(src, destination)
    _write_json(metadata / "wgs_windows.json", n_value_metadata(meta))
    by_key = {(r["contig"], r["haplotype"], int(r["raw_start"]), int(r["actual_length"])): r
              for r in rows if str(r["kept"]).lower() in {"1", "true"}}
    fields = [*N_VALUE_FIELDS, "bl_frac", "coord_coverage", "storage_policy"]
    annotated = 0
    for record in records:
        path = compact_root / record["contig"] / "windows.tsv"
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            names = [name for name in (reader.fieldnames or []) if name not in LEGACY_N_FLAG_FIELDS]
            windows = list(reader)
        for row in windows:
            for field in LEGACY_N_FLAG_FIELDS:
                row.pop(field, None)
            key = (row["contig"], row["haplotype"], int(row["raw_start"]), int(row["actual_length"]))
            if key not in by_key:
                raise ValueError(f"compact/window index 좌표 불일치: {key}")
            row.update({field: by_key[key][field] for field in fields})
            annotated += 1
        part = path.with_name(path.name + ".part")
        with part.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=names + [f for f in fields if f not in names], lineterminator="\n")
            writer.writeheader()
            writer.writerows(windows)
        os.replace(part, path)
    # max_windows는 debug subset일 수 있으므로 전체 candidate 원장도 별도로 보존한다.
    mappings = []
    for record in records:
        for hap in ("hap1", "hap2"):
            raw = record.get(hap + "_chain", "")
            if not raw:
                mode = "unavailable_assembly_direct" if sequence_source == "assembly_direct" else "identity"
                mappings.append({"contig": record["contig"], "haplotype": hap.upper(), "mapping": mode})
                continue
            src = Path(raw)
            destination = metadata / "reference_mapping" / hap.upper() / (record["contig"] + ".chain.gz")
            destination.parent.mkdir(parents=True, exist_ok=True)
            part = destination.with_name(destination.name + ".part")
            with src.open("rb") as inp, part.open("wb") as out:
                with gzip.GzipFile(fileobj=out, mode="wb", mtime=0) as zipped:
                    shutil.copyfileobj(inp, zipped)
            source_hash = _sha256(src)
            check = hashlib.sha256()
            with gzip.open(part, "rb") as handle:
                for block in iter(lambda: handle.read(1 << 20), b""):
                    check.update(block)
            if check.hexdigest() != source_hash:
                raise ValueError(f"최종 chain 복사 검증 실패: {src}")
            os.replace(part, destination)
            mappings.append({"contig": record["contig"], "haplotype": hap.upper(), "mapping": "ucsc_chain",
                             "path": str(destination.relative_to(sample_root)), "source_sha256": source_hash,
                             "stored_sha256": _sha256(destination), "stored_bytes": destination.stat().st_size})
    summary = {"schema": "mfm_window_annotations_v1", "storage_policy": POLICY,
               "n_annotation_schema": N_ANNOTATION_SCHEMA, "n_value_definitions": N_VALUE_DEFINITIONS,
               "terminal_n": "trim continuous chromosome/HAP terminal N; preserve original offsets",
               "internal_n": "preserve positions; exclude windows without any ACGT",
               "blacklist": "annotate only at default max_blacklist_frac=1.0",
               "candidate_windows": len(rows), "stored_windows": annotated,
               "compact_root": str(compact_root.relative_to(sample_root)),
               "candidate_manifest_sha256": _sha256(metadata / "wgs_windows.tsv"),
               "training_selection": {"applied_to_storage": False,
                                      "status": "consumer_defined_from_n_fraction"},
               "dna_target_contract": "consumer must exclude N and padding P target positions; normalize by valid target count",
               "rna_wgs_contract": "use observed alignment/base quality; N/blacklist does not overwrite observed signal",
               "reference_contract": reference_contract, "mappings": mappings}
    _write_json(metadata / "window_policy.json", summary)
    return summary


def validate_window_annotations(sample_root: Path) -> dict:
    """삭제 전과 삭제 후 모두 사용할 수 있는 최종 정책·매핑 무결성 검사."""
    metadata = sample_root / "metadata"
    summary = json.loads((metadata / "window_policy.json").read_text())
    if summary.get("storage_policy") != POLICY:
        raise ValueError("window 저장 정책 식별자 불일치")
    annotation_schema = summary.get("n_annotation_schema")
    if annotation_schema not in (None, "mfm_window_n_annotations_v2", N_ANNOTATION_SCHEMA):
        raise ValueError("지원하지 않는 window N 주석 schema")
    # 과거 최종 패키지는 읽기 전용 호환을 유지함. 새 schema는 원값만 저장함.
    if _sha256(metadata / "wgs_windows.tsv") != summary["candidate_manifest_sha256"]:
        raise ValueError("최종 후보 원장 SHA256 불일치")
    with (metadata / "wgs_windows.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != summary["candidate_windows"]:
        raise ValueError("최종 후보 창 개수 불일치")
    for row in rows:
        validate_n_values(row)
        if annotation_schema == N_ANNOTATION_SCHEMA and set(LEGACY_N_FLAG_FIELDS).intersection(row):
            raise ValueError("원값 schema에 파생 N flag가 남아 있음")
        if str(row["kept"]).lower() in {"1", "true"} and int(row["n_nonacgt"]) >= int(row["actual_length"]):
            raise ValueError("모든 염기가 미상인 창이 최종 입력에 포함됨")
    by_key = {(r["contig"], r["haplotype"], int(r["raw_start"]), int(r["actual_length"])): r
              for r in rows if str(r["kept"]).lower() in {"1", "true"}}
    observed = 0
    compact_root = sample_root / summary["compact_root"]
    if not compact_root.resolve().is_relative_to(sample_root.resolve()):
        raise ValueError("compact 주석 경로가 sample 밖을 가리킴")
    for path in compact_root.glob("*/windows.tsv"):
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if annotation_schema == N_ANNOTATION_SCHEMA and set(LEGACY_N_FLAG_FIELDS).intersection(row):
                    raise ValueError("compact 원값 schema에 파생 N flag가 남아 있음")
                key = (row["contig"], row["haplotype"], int(row["raw_start"]), int(row["actual_length"]))
                if key not in by_key:
                    raise ValueError("compact 주석의 후보 좌표 불일치")
                for field in (*N_VALUE_FIELDS, "bl_frac", "coord_coverage", "storage_policy"):
                    if row.get(field) != by_key[key].get(field):
                        raise ValueError(f"compact 주석 불일치: {field}")
                observed += 1
    if observed != summary["stored_windows"]:
        raise ValueError("compact 주석의 저장 창 개수 불일치")
    for mapping in summary["mappings"]:
        if mapping["mapping"] != "ucsc_chain":
            continue
        path = sample_root / mapping["path"]
        if not path.resolve().is_relative_to(sample_root.resolve()):
            raise ValueError("최종 mapping 경로가 sample 밖을 가리킴")
        if _sha256(path) != mapping["stored_sha256"]:
            raise ValueError("최종 mapping 파일 SHA256 불일치")
        digest = hashlib.sha256()
        with gzip.open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        if digest.hexdigest() != mapping["source_sha256"]:
            raise ValueError("최종 mapping 원본 SHA256 불일치")
    return {"status": "pass", "candidate_windows": len(rows), "stored_windows": summary["stored_windows"],
            "mapping_records": len(summary["mappings"])}
