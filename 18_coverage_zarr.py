#!/usr/bin/env python3
"""정렬 BAM → 전장 및 WGS-window 1 bp coverage Zarr.

GTEx `bam2coverage.py`와 같이 MAPQ255 primary·QC-pass alignment를 사용하고
duplicate는 포함한다. ASE fragment count의 duplicate 제외 정책과는 목적이 다르다.

    python 18_coverage_zarr.py --bam <bam> --fai <target.fa.fai> --out <coverage.zarr> \
        --windows-tsv <windows.tsv> --window-out <coverage_windows.zarr> \
        --contigs "chr1 chr2 ..." --sample S --hap HAP1

★ read 선택: MAPQ255 primary QC-pass이고 secondary/supplementary를 제외한다.
  `samtools depth ... -g 1024`로 기본 DUP 제외를 해제한다. splice N 제외,
  deletion 포함, mate 겹침 2회 계수이다.
★ 정규화: **원본 depth는 정규화하지 않고 보존한다.** 정규화 상수만 zarr attrs 에 적고 변환은 dataloader 가 한다 —
  지수·clip 임계를 바꿔 재실험할 수 있어야 하고, ASE 처럼 비율을 보는 분석은 raw 가 필요하다.
      m = nonzero_mean = mean(x | x>0)
      y = soft_clip((x/m)^power, clip)        power=0.75, clip=10  (AlphaGenome)
  `nonzero_mean` 은 두 범위로 저장한다(genome 전체 / 나중에 kept window 로 재계산 가능).
"""
from __future__ import annotations
import argparse, json, os, sys, tempfile
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc,blosc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rna_profile  # noqa: E402
from rna_coverage import collect_coverage  # noqa: E402
from rna_resources import file_signature  # noqa: E402
from rna_window_coverage import write_window_coverage  # noqa: E402

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--windows-tsv", type=Path)
    ap.add_argument("--window-out", type=Path)
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--chunk", type=int, default=4194304)
    ap.add_argument("--power", type=float, default=0.75)
    ap.add_argument("--clip", type=float, default=10.0)
    ap.add_argument("--sample", default="")
    ap.add_argument("--hap", default="")
    ap.add_argument("--samtools", default=os.environ.get("SAMTOOLS", "samtools"))
    ap.add_argument("--jobs", default=os.environ.get("RNA_COVERAGE_JOBS", "auto"))
    a = ap.parse_args()
    if bool(a.windows_tsv) != bool(a.window_out):
        ap.error("--windows-tsv와 --window-out은 함께 지정해야 합니다")

    lengths = {}
    for line in open(a.fai):
        f = line.split("\t"); lengths[f[0]] = int(f[1])
    contigs = [c for c in a.contigs.split() if c in lengths]
    if len(set(contigs)) != len(contigs):
        raise ValueError("병렬 coverage 입력에 contig가 중복됐습니다")
    source_signature = file_signature(a.bam)

    blosc.set_nthreads(1)
    comp = Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE)
    g = zarr.open(a.out, mode="w")
    for c in contigs:
        L = lengths[c]
        g.create_dataset(c, shape=(L,), dtype=np.uint32, chunks=(min(a.chunk, L),),
                         compressor=comp, overwrite=True)
    with tempfile.TemporaryDirectory(prefix='rna_buffers_') as buffers:
        try:
            nz_sum, nz_n, total, depth_sha256 = collect_coverage(
                contigs, lengths, a.bam, a.samtools, a.out, a.jobs, snapshot_dir=buffers if a.windows_tsv else None
            )
        except RuntimeError as exc:
            print(f"[FAIL] {exc}", file=sys.stderr, flush=True)
            return 1
        if file_signature(a.bam) != source_signature:
            raise RuntimeError("coverage 계산 중 source BAM이 변경됐습니다")

        m = nz_sum / nz_n if nz_n else 0.0
        if total <= 0 or nz_n <= 0:
            print(f"[FAIL] {a.sample} {a.hap}: coverage가 0입니다", file=sys.stderr)
            return 1
        g.attrs.update({
            "encoding": ("unscaled read depth, uint32, 1 bp resolution "
                         "(samtools depth -a -J: deletion 포함, splice N 미포함, mate 겹침 2회 계수)"),
            "mate_overlap_policy": "counted twice (samtools depth -s not used)",
            "deletion_policy": "counted as covered (samtools depth -J)",
            "sample": a.sample, "hap": a.hap, "source": "STAR alignment to personal haplotype",
            "coordinate_space": "individual (personal haplotype)",
            "nonzero_mean": m, "nonzero_mean_scope": "genome",
            "norm_power": a.power, "norm_clip_at": a.clip,
            "norm_recipe": "AlphaGenome: y = soft_clip((x / nonzero_mean) ** power, clip_at); raw signal 은 그대로 보존한다",
            "total_depth": total, "nonzero_positions": nz_n,
            "depth_sha256_encoding": "per-contig raw little-endian uint32 C-order",
            "depth_sha256_by_contig": depth_sha256,
            "min_mapq_applied": 255,
            "source_bam_signature": source_signature,
            "source_fai_signature": file_signature(a.fai),
            **rna_profile.coverage_metadata(),
        })
        if os.environ.get('MFM_COVERAGE_BACKEND','cpu')=='cupy':
            g.attrs['compute_execution']={'coverage_backend':'cupy','decode_backend':'HTSlib/native CPU; GPU child exits before Zarr compression','storage_backend':'Blosc-Zstd5/CPU'}
        print(f"[cov] {a.sample} {a.hap}: 총 depth {total:,} / >0 위치 {nz_n:,} / nonzero_mean {m:.4f}")
        if a.windows_tsv:
            window_stats = write_window_coverage(
                source_path=Path(a.out), windows_path=a.windows_tsv, output_path=a.window_out,
                sample=a.sample, hap=a.hap, lengths=lengths, contigs=contigs,
                chunk=a.chunk, power=a.power, clip=a.clip, source_buffers=Path(buffers),
            )
            print("[cov:windows] " + json.dumps(window_stats, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
