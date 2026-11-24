#!/usr/bin/env python3
"""[구번호] RNA 산출물을 `metadata/sample.json` 에 반영한다.

현재 v2 대응 단계는 Stage 24 (`24_stamp_rna_metadata.py`)다. 기존 metadata 진입점과
stamp helper만 공유하며, Stage 24의 QC 경로 기록 범위가 더 넓다.

    python 19_stamp_rna_metadata.py <final_sample_root> [<final_sample_root> ...] [--dry-run]

★ 왜 별도 단계인가
  `sample.json` 은 WGS 통합 전처리(`preprocess-full-shared-hap`)가 쓴다. 그런데 RNA 파이프라인은
  그보다 **뒤에 도는 경우가 있어서**, 그 시점에는 `rna_personal/` 이 아직 없다. 그래서 WGS 쪽에서
  아무리 잘 판정해도 `has_rna` 가 false 로 남을 수 있다. RNA 가 끝난 시점에 한 번 더 찍어야 한다.

★ 무엇을 고치나
  - `has_rna` : bigwig 또는 rna_personal 이 있으면 true
  - `rna_personal` : `{rna_personal_root, phase_mode, coverage_zarr{HAP1,HAP2}, allelic_counts}`
    로더가 경로를 추측하지 않게 실제 경로를 실어 보낸다. 신규 경로에서 RNA track 은
    `tracks/` 가 아니라 `rna_personal/<phase>/HAP{1,2}/coverage.zarr` 에 있다.

  배경: `post-training/journal/2609/post_training_data_issue_260903.md` §5·§6
  (신규 10명 전부 `has_rna: false` 로 잘못 박혀 있었다 — 이 필드로 필터링하면 RNA 가 통째로 누락된다)

기존 키는 건드리지 않고 위 두 개만 갱신한다. 여러 번 돌려도 결과가 같다.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


from rna_metadata_common import detect_rna_personal as _detect_rna_personal, stamp_sample

def detect_rna_personal(sample_root: Path) -> dict:
    return _detect_rna_personal(sample_root, include_variant_qc=False)

def stamp(sample_root: Path, dry_run: bool) -> str:
    return stamp_sample(sample_root, dry_run, detect_rna_personal)

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("sample_roots", nargs="+", help="final/.../<cohort>/<sample> 경로")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    rc = 0
    for raw in a.sample_roots:
        p = Path(raw)
        if not p.is_dir():
            print(f"FAIL  {raw}: 디렉터리 없음", file=sys.stderr); rc = 1; continue
        print(stamp(p, a.dry_run))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
