#!/usr/bin/env python3
"""[구번호] 정렬 BAM → 1 bp 커버리지 zarr (개인 좌표) + AlphaGenome 정규화 상수.

현재 v2 runner의 대응은 Stage 18 (`18_coverage_zarr.py`)이다. 이 파일의 filtering
설명과 실행 결과는 구 구현 기준이므로 새 전처리의 정책 설명으로 사용하지 않는다.

    python 15_coverage_zarr.py --bam <bam> --fai <target.fa.fai> --out <coverage.zarr> \
        --contigs "chr1 chr2 ..." --sample S --hap HAP1

★ 규칙 (instruction/preprocessing_rule.md "RNA 타깃 정규화 규칙")
  **원본 depth는 정규화하지 않고 보존한다.** STAR 뒤 Picard MarkDuplicates가 표기한
  duplicate는 `samtools depth` 기본 필터로 제외한다. 필요 시에만 `--include-duplicates`로
  명시적으로 포함한다. 정규화 상수만 zarr attrs 에 적고 변환은 dataloader 가 한다 —
  지수·clip 임계를 바꿔 재실험할 수 있어야 하고, ASE 처럼 비율을 보는 분석은 raw 가 필요하다.
      m = nonzero_mean = mean(x | x>0)
      y = soft_clip((x/m)^power, clip)        power=0.75, clip=10  (AlphaGenome)
  `nonzero_mean` 은 두 범위로 저장한다(genome 전체 / 나중에 kept window 로 재계산 가능).
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys

import numpy as np
import zarr
from numcodecs import Blosc

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--chunk", type=int, default=4194304)
    ap.add_argument("--power", type=float, default=0.75)
    ap.add_argument("--clip", type=float, default=10.0)
    ap.add_argument("--sample", default="")
    ap.add_argument("--hap", default="")
    ap.add_argument("--min-mapq", type=int, default=0,
                    help="0 = 필터 없음. 파이프라인 원칙: 필터는 소비 지점에서 (기본 0 을 바꾸지 말 것)")
    ap.add_argument("--include-duplicates", action="store_true",
                    help="Picard DUP flag를 포함한다. 기본값은 samtools depth의 DUP 제외 정책이다")
    ap.add_argument("--samtools", default=os.environ.get("SAMTOOLS", "samtools"))
    a = ap.parse_args()

    lengths = {}
    for line in open(a.fai):
        f = line.split("\t"); lengths[f[0]] = int(f[1])
    contigs = [c for c in a.contigs.split() if c in lengths]

    comp = Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE)
    g = zarr.open(a.out, mode="w")
    nz_sum = 0; nz_n = 0; total = 0
    for c in contigs:
        L = lengths[c]
        cov = np.zeros(L, dtype=np.uint32)
        # samtools depth -a -J: 모든 위치(-a), splice(N) 구간은 커버리지에서 제외(-J)
        cmd = [a.samtools, "depth", "-a", "-J", "-r", c]
        if a.include_duplicates:
            cmd += ["-g", "DUP"]
        if a.min_mapq:
            cmd += ["-Q", str(a.min_mapq)]
        cmd.append(a.bam)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.split("\t")
            i = int(f[1]) - 1
            if 0 <= i < L:
                cov[i] = int(f[2])
        p.wait()
        if p.returncode != 0:
            print(f"[FAIL] samtools depth {c} exit {p.returncode}", file=sys.stderr); return 1
        g.create_dataset(c, data=cov, chunks=(min(a.chunk, L),), compressor=comp, overwrite=True)
        nzm = cov > 0
        nz_sum += int(cov[nzm].sum()); nz_n += int(nzm.sum()); total += int(cov.sum())
        print(f"  {c}: len {L:,} / 합 {int(cov.sum()):,} / >0 {int(nzm.sum()):,}")

    m = nz_sum / nz_n if nz_n else 0.0
    if total <= 0 or nz_n <= 0:
        print(f"[FAIL] {a.sample} {a.hap}: coverage가 0입니다", file=sys.stderr)
        return 1
    g.attrs.update({
        "encoding": "unscaled read depth, uint32, 1 bp resolution (splice N 구간 제외: samtools depth -J)",
        "sample": a.sample, "hap": a.hap, "source": "STAR alignment to personal haplotype",
        "coordinate_space": "individual (personal haplotype)",
        "nonzero_mean": m, "nonzero_mean_scope": "genome",
        "norm_power": a.power, "norm_clip_at": a.clip,
        "norm_recipe": "AlphaGenome: y = soft_clip((x / nonzero_mean) ** power, clip_at); raw signal 은 그대로 보존한다",
        "total_depth": total, "nonzero_positions": nz_n,
        "min_mapq_applied": a.min_mapq,
        "duplicate_policy": "included (--include-duplicates)" if a.include_duplicates else "excluded (samtools depth default DUP filter)",
    })
    print(f"[cov] {a.sample} {a.hap}: 총 depth {total:,} / >0 위치 {nz_n:,} / nonzero_mean {m:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
