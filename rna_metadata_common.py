#!/usr/bin/env python3
"""RNA 산출물을 sample metadata에 기록하는 공통 helper.

Stage 19/23 구번호와 현재 Stage 24는 같은 stamp 동작을 쓴다. Stage 24만
hap-variant QC 경로를 추가 기록하므로 detector의 선택 인자로 차이를 보존한다.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Callable

_BASE_OUTPUTS = {
    "allelic_counts_meta": "allelic_counts.meta.json",
    "allelic_counts_split_half": "allelic_counts.split_half.tsv",
    "allelic_counts_both_haps": "allelic_counts.both_haps.tsv",
    "ase_qc_tsv": "ase_qc.tsv",
    "ase_qc_json": "ase_qc.json",
    "ase_training_eligibility": "ase_training_eligibility.tsv",
}
_VARIANT_OUTPUTS = {
    "hap_variants": "hap_variants.tsv",
    "hap_variants_qc": "hap_variants_qc.json",
}


def detect_rna_personal(sample_root: Path, *, include_variant_qc: bool = False) -> dict[str, Any]:
    """`rna_personal/<phase>/`의 실제 결과 경로를 찾아 반환한다."""
    root = sample_root / "rna_personal"
    if not root.is_dir():
        return {}
    optional_outputs = dict(_BASE_OUTPUTS)
    if include_variant_qc:
        optional_outputs = {**_VARIANT_OUTPUTS, **optional_outputs}
    for phase_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        haps = {
            hap: str(phase_dir / hap / "coverage.zarr")
            for hap in ("HAP1", "HAP2")
            if (phase_dir / hap / "coverage.zarr").is_dir()
        }
        if not haps:
            continue
        out: dict[str, Any] = {
            "rna_personal_root": str(phase_dir),
            "phase_mode": phase_dir.name,
            "coverage_zarr": haps,
        }
        allelic = phase_dir / "allelic_counts.tsv"
        if allelic.is_file():
            out["allelic_counts"] = str(allelic)
        for key, filename in optional_outputs.items():
            path = phase_dir / filename
            if path.is_file():
                out[key] = str(path)
        return out
    return {}


def stamp_sample(
    sample_root: Path,
    dry_run: bool,
    detector: Callable[[Path], dict[str, Any]],
) -> str:
    """기존 metadata를 보존하며 RNA 상태만 원자적으로 갱신한다."""
    meta = sample_root / "metadata" / "sample.json"
    if not meta.is_file():
        return f"SKIP  {sample_root.name}: metadata/sample.json 없음"
    data = json.loads(meta.read_text())
    rna = detector(sample_root)
    has_rna = bool(data.get("rna_bigwigs")) or bool(rna)
    before = (data.get("has_rna"), bool(data.get("rna_personal")))
    if data.get("has_rna") == has_rna and data.get("rna_personal") == rna:
        return f"OK    {sample_root.name}: 이미 최신 (has_rna={has_rna})"
    data["has_rna"] = has_rna
    data["rna_personal"] = rna
    if dry_run:
        return f"DRY   {sample_root.name}: has_rna {before[0]} → {has_rna}, rna_personal 기록"
    tmp = meta.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    os.replace(tmp, meta)
    return f"WROTE {sample_root.name}: has_rna {before[0]} → {has_rna}, phase={rna.get('phase_mode','-')}"
