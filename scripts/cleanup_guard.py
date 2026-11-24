#!/usr/bin/env python3
"""중간 산출물 삭제 직전에 최종 검증과 RNA ASE provenance를 다시 확인한다."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


FINAL_SCHEMA = "mfm_final_package_validation_v1"
FINAL_REVISION = "2026-10-03_integrity_v2"
HAP_QC_SCHEMA = "mfm_hap_variants_qc_v1"
COUNT_SCHEMA = "mfm_allelic_counts_meta_v5"
COUNT_REVISION = "2026-10-02_validated_hap_variant_evidence"
ASE_SCHEMA = "mfm_ase_qc_v3"
ASE_REVISION = "2026-10-02_validated_hap_variant_evidence_cnv_hap_consistency"

RNA_REQUIRED_FILES = (
    "hap_variants.tsv",
    "hap_variants_qc.json",
    "allelic_counts.tsv",
    "allelic_counts.meta.json",
    "allelic_counts.split_half.tsv",
    "allelic_counts.both_haps.tsv",
    "allelic_counts.mmap_adjusted.tsv",
    "ase_qc.json",
    "ase_training_eligibility.tsv",
    "rna_realign_summary.json",
)


def _file(path: Path) -> Path:
    if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"삭제 전 필수 파일 없음/비정상: {path}")
    return path


def _directory(path: Path) -> Path:
    if path.is_symlink() or not path.is_dir():
        raise ValueError(f"삭제 전 필수 디렉터리 없음/비정상: {path}")
    return path


def _json(path: Path) -> dict:
    try:
        value = json.loads(_file(path).read_text())
    except Exception as exc:
        raise ValueError(f"삭제 전 JSON 판독 실패: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"삭제 전 JSON object 형식 아님: {path}")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


def validate_predelete(root: str | Path, sample: str, mode: str) -> dict:
    """현재 파일을 직접 검사한다. 예전 PASS JSON만으로 삭제를 허가하지 않는다."""
    root = Path(root).absolute()
    if mode not in {"wgs", "rna"}:
        raise ValueError(f"삭제 전 mode 오류: {mode}")
    if root.is_symlink() or not root.is_dir() or root.name != sample:
        raise ValueError(f"삭제 전 sample root 불일치: root={root} sample={sample}")

    report_path = root / "qc" / "final_package_validation.json"
    report = _json(report_path)
    if (report.get("schema") != FINAL_SCHEMA or report.get("contract_revision") != FINAL_REVISION
            or report.get("sample") != sample
            or report.get("mode") != mode):
        raise ValueError("삭제 전 최종 검증 schema/revision/sample/mode 불일치")
    if report.get("status") != "pass" or report.get("failed_checks"):
        raise ValueError("삭제 전 최종 검증 PASS가 필요함")

    personal_manifest = _file(root / "individual_fasta" / "permanent_assets.json")
    if report.get("personal_assets_manifest_sha256") != _sha256(personal_manifest):
        raise ValueError("삭제 전 영구 HAP 자산 manifest SHA-256 불일치")

    result = {
        "schema": "mfm_cleanup_predelete_guard_v1",
        "sample": sample,
        "mode": mode,
        "validation_report": str(report_path),
        "personal_assets_manifest": str(personal_manifest),
    }
    if mode == "wgs":
        return result

    rna_base = _directory(root / "rna_personal")
    labels = sorted(p for p in rna_base.iterdir() if p.is_dir() and not p.is_symlink())
    if len(labels) != 1:
        raise ValueError(f"삭제 전 RNA evidence 디렉터리는 하나여야 함: {[p.name for p in labels]}")
    rna = labels[0]
    if report.get('rna_payload','full')=='coverage_only':
        from rna_signal import validate_signal
        validate_signal(root,sample,rna.name)
        if report.get('rna_signal_manifest_sha256')!=_sha256(rna/'rna_signal_manifest.json'):
            raise ValueError('삭제 전 RNA signal manifest SHA 불일치')
        newest=max(p.stat().st_mtime_ns for p in rna.rglob('*') if p.is_file())
        if report_path.stat().st_mtime_ns<newest:raise ValueError('삭제 전 RNA signal 검증 보고서가 오래됨')
        result.update(rna_payload='coverage_only',rna_label=rna.name,validation_report_fresh=True)
        return result
    files = {name: _file(rna / name) for name in RNA_REQUIRED_FILES}
    for hap in ("HAP1", "HAP2"):
        _directory(rna / hap / "coverage.zarr")
        _directory(rna / hap / "coverage_windows.zarr")

    hap_sha = _sha256(files["hap_variants.tsv"])
    hap_qc = _json(files["hap_variants_qc.json"])
    if hap_qc.get("schema") != HAP_QC_SCHEMA or hap_qc.get("output_sha256") != hap_sha:
        raise ValueError("삭제 전 HAP variant QC schema/SHA-256 불일치")

    count_meta = _json(files["allelic_counts.meta.json"])
    count_paths = {
        "counts": files["allelic_counts.tsv"],
        "split_half": files["allelic_counts.split_half.tsv"],
        "both_haps": files["allelic_counts.both_haps.tsv"],
        "mmap_adjusted": files["allelic_counts.mmap_adjusted.tsv"],
    }
    observed_count_hashes = {name: _sha256(path) for name, path in count_paths.items()}
    if (count_meta.get("schema") != COUNT_SCHEMA
            or count_meta.get("script_revision") != COUNT_REVISION
            or count_meta.get("hap_variants_sha256") != hap_sha
            or count_meta.get("output_sha256") != observed_count_hashes):
        raise ValueError("삭제 전 allelic count schema/revision/HAP/output SHA-256 불일치")

    ase_qc = _json(files["ase_qc.json"])
    sample_qc = ase_qc.get("sample_qc", {})
    inputs = ase_qc.get("inputs", {})
    input_hash_ok = (
        inputs.get("counts_sha256") == observed_count_hashes["counts"]
        and inputs.get("split_half_sha256") == observed_count_hashes["split_half"]
        and inputs.get("both_haps_sha256") == observed_count_hashes["both_haps"]
        and inputs.get("mmap_adjusted_sha256") == observed_count_hashes["mmap_adjusted"]
        and inputs.get("meta_sha256") == _sha256(files["allelic_counts.meta.json"])
    )
    cnv_path_raw = inputs.get("cnv_rd")
    if cnv_path_raw:
        cnv_path = _file(Path(cnv_path_raw))
        cnv_meta_path = _file(Path(str(cnv_path) + ".meta.json"))
        cnv_meta = _json(cnv_meta_path)
        input_hash_ok = (
            input_hash_ok
            and inputs.get("cnv_rd_sha256") == _sha256(cnv_path)
            and inputs.get("cnv_meta_sha256") == _sha256(cnv_meta_path)
            and cnv_meta.get("output_sha256") == _sha256(cnv_path)
        )
    if (ase_qc.get("schema") != ASE_SCHEMA
            or ase_qc.get("script_revision") != ASE_REVISION
            or inputs.get("hap_variants_sha256") != hap_sha
            or not input_hash_ok
            or sample_qc.get("provenance_ok") is not True
            or sample_qc.get("duplicates_excluded") is not True):
        raise ValueError("삭제 전 ASE schema/revision/input SHA-256/provenance 불일치")

    newest_evidence = max(path.stat().st_mtime_ns for path in rna.rglob("*") if path.is_file())
    if report_path.stat().st_mtime_ns < newest_evidence:
        raise ValueError("삭제 전 최종 검증 보고서가 RNA evidence보다 오래됨; Stage 검증 재실행 필요")

    result.update({
        "rna_label": rna.name,
        "hap_variants_sha256": hap_sha,
        "rna_files_checked": len(files),
        "validation_report_fresh": True,
    })
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="cleanup 삭제 전 fail-closed guard")
    parser.add_argument("--root", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mode", required=True, choices=("wgs", "rna"))
    args = parser.parse_args()
    print(json.dumps(validate_predelete(args.root, args.sample, args.mode), ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
