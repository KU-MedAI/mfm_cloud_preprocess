#!/usr/bin/env python3
"""[구번호] WGS 깊이로 유전자별 copy-number 비율(rd)을 재고 CNV 의심 구간을 표시한다.

현재 v2 대응은 Stage 22 (`22_cnv_filter.py`)이며 이 파일은 과거 단일 작업 방식으로 보존한다.

    python 17_cnv_filter.py --exons <exons.bygene.tsv> --bam <personal.HAPn.cram> \
        --out <gene_rd.tsv> [--reference <HAPn.fa>] [--samtools samtools]

★ 왜 필요한가 (AlleleSeq2 대조)
  CNV 영역의 유전자는 **구조적으로** allele-specific 하게 보인다. 한쪽 haplotype 의 사본 수가
  많으면 그쪽 리드가 많은 것이 당연한데, 전사 조절 차이로 읽히면 위양성이다.
  1000G 처럼 다양한 집단이 섞인 코호트에서는 확실한 위양성 소스다.

  AlleleSeq2 는 이를 `makePersonalGenome.mk` 의 `hetSNVs_rd.tab` 타깃에서 처리한다:
  WGS BAM 으로 hetSNV ±1kb 의 read depth 를 재고(`calculate_rd.py`),
  `filter_cnv_sites.py` 가 **0.5 ≤ rd ≤ 1.5 밖이면 그 site 를 버린다.**

  v2 는 site 가 아니라 gene 단위로 세므로 같은 판정을 gene 단위로 옮긴다.
  깊이는 개인 좌표 WGS 정렬에서 잰다 — exon 좌표와 같은 좌표계여야 하므로
  GRCh38 BAM 을 쓰면 안 된다.

★ rd 의 정의
  rd = (그 유전자 exon 합집합의 평균 depth) / (전체 유전자 평균 depth 의 중앙값)
  중앙값으로 나누므로 sample 의 전체 커버리지 깊이에는 영향받지 않는다.

★ duplicate
  `samtools bedcov` 기본 필터를 그대로 쓴다(UNMAP/SECONDARY/QCFAIL/DUP 제외).
  현재 v2에서는 Stage 18 coverage가 duplicate를 포함하고 Stage 21 ASE count는 제외한다.
  이 CNV depth는 samtools bedcov 기본 필터에 따라 DUP를 제외한다.
"""
from __future__ import annotations

import argparse
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path
from rna_interval_common import merge


def load_exons(path: Path) -> dict[str, dict[str, list[tuple[int, int]]]]:
    """gene_id → contig → [(start, end)]. 현재 Stage 14 annotation lift의 4열 형식."""
    per_gene: dict[str, dict[str, list[tuple[int, int]]]] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                continue
            gid, contig, start, end = parts
            per_gene.setdefault(gid, {}).setdefault(contig, []).append((int(start), int(end)))
    return per_gene


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exons", required=True, action="append",
                    help="Stage 14 lift 산출물 exons.bygene.tsv. --bam 과 짝을 맞춰 "
                         "HAP 마다 하나씩 반복 지정한다")
    ap.add_argument("--bam", required=True, action="append",
                    help="Stage 10의 개인 좌표 WGS BAM/CRAM. --exons 와 같은 순서")
    ap.add_argument("--out", required=True)
    ap.add_argument("--reference", default=[], action="append",
                    help="CRAM 디코드용 개인 FASTA. --bam 과 같은 순서(생략 가능)")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--min-rd", type=float, default=0.5, help="AlleleSeq2 filter_cnv_sites.py 와 동일")
    ap.add_argument("--max-rd", type=float, default=1.5)
    ap.add_argument("--min-gene-length", type=int, default=100,
                    help="이보다 짧은 유전자는 깊이 추정이 불안정하므로 rd 를 내지 않는다")
    args = ap.parse_args()

    if len(args.exons) != len(args.bam):
        raise ValueError(f"--exons {len(args.exons)}개와 --bam {len(args.bam)}개의 수가 다릅니다")
    refs = list(args.reference or [])
    if refs and len(refs) != len(args.bam):
        raise ValueError(f"--reference 를 줄 때는 --bam 과 같은 수여야 합니다")

    per_hap_rd: list[dict[str, float]] = []
    per_hap_median: list[float] = []
    gene_length_any: dict[str, int] = {}

    for idx, (exons_arg, bam_arg) in enumerate(zip(args.exons, args.bam)):
        exons_path, bam_path = Path(exons_arg), Path(bam_arg)
        if not exons_path.is_file():
            raise FileNotFoundError(f"exon 파일 없음: {exons_path}")
        if not bam_path.is_file():
            raise FileNotFoundError(f"WGS 정렬 없음: {bam_path}")
        per_gene = load_exons(exons_path)
        if not per_gene:
            raise ValueError(f"exon 이 비어 있습니다: {exons_path}")

        # BED 는 gene 별 병합 구간. 4열에 gene_id 를 넣어 bedcov 출력에서 되찾는다.
        rows: list[tuple[str, int, int, str]] = []
        gene_length: dict[str, int] = {}
        for gid, by_contig in per_gene.items():
            for contig, ivs in by_contig.items():
                for start, end in merge(ivs):
                    if end <= start:
                        continue
                    rows.append((contig, start, end, gid))
                    gene_length[gid] = gene_length.get(gid, 0) + (end - start)
        rows.sort(key=lambda r: (r[0], r[1], r[2]))
        for gid, length in gene_length.items():
            gene_length_any[gid] = max(gene_length_any.get(gid, 0), length)

        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "genes.bed"
            with bed.open("w") as fh:
                for contig, start, end, gid in rows:
                    fh.write(f"{contig}\t{start}\t{end}\t{gid}\n")
            cmd = [args.samtools, "bedcov"]
            if refs and refs[idx]:
                cmd += ["--reference", refs[idx]]
            cmd += [str(bed), str(bam_path)]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                print(f"[FAIL] samtools bedcov exit {proc.returncode}\n{proc.stderr[-2000:]}",
                      file=sys.stderr)
                return 1
            depth_sum: dict[str, float] = {}
            for line in proc.stdout.splitlines():
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                gid, covered = parts[3], parts[4]
                depth_sum[gid] = depth_sum.get(gid, 0.0) + float(covered)

        mean_depth: dict[str, float] = {}
        for gid, length in gene_length.items():
            if length >= args.min_gene_length:
                mean_depth[gid] = depth_sum.get(gid, 0.0) / length
        if not mean_depth:
            raise ValueError(f"깊이를 잴 수 있는 유전자가 없습니다: {bam_path}")
        median_depth = statistics.median(mean_depth.values())
        if median_depth <= 0:
            raise ValueError(
                f"유전자 평균 depth 의 중앙값이 0 입니다 — WGS 정렬과 exon 좌표계가 다를 수 있습니다: {bam_path}")
        per_hap_rd.append({g: d / median_depth for g, d in mean_depth.items()})
        per_hap_median.append(median_depth)
        print(f"[cnv] hap{idx + 1}: 유전자 {len(gene_length):,} / rd 산출 {len(mean_depth):,} / "
              f"중앙 depth {median_depth:.2f}x")

    n_hap = len(per_hap_rd)
    n_pass = n_fail = n_na = 0
    out_path = Path(args.out)
    with out_path.open("w") as fh:
        hap_cols = "".join(f"\trd_hap{i + 1}" for i in range(n_hap))
        fh.write(f"gene_id\texon_bp{hap_cols}\trd\tpass_cnv\n")
        for gid in sorted(gene_length_any):
            vals = [hap.get(gid) for hap in per_hap_rd]
            cells = "".join(f"\t{v:.4f}" if v is not None else "\tNA" for v in vals)
            present = [v for v in vals if v is not None]
            if not present:
                # 판정할 깊이가 없다(너무 짧거나 한쪽 hap 에만 있음) → 통과로 두고 NA 로 남긴다.
                n_na += 1
                fh.write(f"{gid}\t{gene_length_any[gid]}{cells}\tNA\t1\n")
                continue
            rd = statistics.fmean(present)
            # 한쪽 hap 이라도 구간을 벗어나면 제외한다(보수적).
            ok = all(args.min_rd <= v <= args.max_rd for v in present)
            n_pass += ok
            n_fail += not ok
            fh.write(f"{gid}\t{gene_length_any[gid]}{cells}\t{rd:.4f}\t{1 if ok else 0}\n")

    print(f"[cnv] 통과 {n_pass:,} / CNV 의심 제외 {n_fail:,} / 판정보류 {n_na:,} "
          f"(구간 {args.min_rd}–{args.max_rd}, hap {n_hap}개 모두 통과 요구) → {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
