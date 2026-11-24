#!/usr/bin/env python3
"""[구번호] 두 haplotype 정렬을 **fragment 별 AS 비교**로 allele count를 만든다.

현재 v2 Stage 21은 `21_allelic_counts.py`이며 Stage 20의 검증된 HAP variant 증거를
요구한다. 이 파일은 variant evidence 적용 전 계산의 재현·비교용으로 보존한다.

    python 20_allelic_counts.py --bam1 <HAP1.bam> --bam2 <HAP2.bam> \
        --exons1 <exons.hap1.tsv> --exons2 <exons.hap2.tsv> --out <allelic_counts.tsv>

★ 왜 이게 필요한가
  H1 에 정렬한 커버리지는 **총 발현**이다(H2 리드도 H1 에 붙는다). haplotype 특이 신호는
  같은 조각이 두 haplotype 중 **어느 쪽에 더 잘 붙는지**로만 나온다:
      AS(H1) > AS(H2)  → hap1
      AS(H1) < AS(H2)  → hap2
      AS(H1) == AS(H2) → **무정보(버린다)** — het 를 지나지 않은 조각이다
  실측(HG00096): het 평균 간격 1,067 bp / exon 안 1,441 bp ⇒ 75 bp 리드의 5.1%,
  ~150 bp fragment 의 9.9% 만 배정 가능하다. 이건 방법의 한계가 아니라 het 밀도의 물리적 상한이고
  phASER 도 같은 제약을 받는다. **재정렬로 개선되는 것은 배정량이 아니라 reference mapping bias 다.**

★ 2026-09-04 수정 — read 단위 → **fragment 단위**
  이전 구현은 키를 `(query_name, read1|read2)` 로 잡아 **한 fragment 의 두 mate 를 독립된 증거
  2개로 셌다.** 두 mate 는 같은 분자라 haplotype 이 반드시 같으므로, 표본수만 2배가 되고
  정보량은 그대로다. 결과적으로 `totalCount` 가 부풀고 이항검정이 과도하게 유의해졌다 —
  실측으로 `is_ase`(FDR<0.05) 가 10명 평균 **0.45** 였다(참조축 진리값 0.067, LCL 정상 범위 5–10%).
  근거: `post-training/journal/2609/ase_truth_spread_260903.md` §9.4

  이제 fragment 하나가 한 표다. 점수는 **양쪽에 공통으로 정렬된 mate 에서만** 합산해 비교한다 —
  한쪽에만 있는 mate 를 합에 넣으면 mate 수 차이가 그대로 점수 차이가 되어 "정렬 성공/실패"가
  대립 증거로 둔갑한다. 유전자 배정은 이긴 haplotype 쪽 mate 들의 blocks **합집합**으로 하고
  유전자당 +1 이다(이전에는 mate 마다 +1 이라 한 유전자가 +2 를 받을 수 있었다).

★ 한쪽만 정렬된 fragment (only1/only2)
  AS 비교가 아니라 정렬 성공/실패에서 나온 신호라 대립 증거로 보기 약하다. 전체 배정의 10–19%.
  다만 이를 빼는 것은 별개의 판단이므로 **기본 동작은 유지**하고 `--require-both-haps` 로만
  제외할 수 있게 한다. 한 번에 한 변수만 바꿔야 위 수정의 효과를 검증할 수 있다.
  ⚠ AlleleSeq2 는 이 구분 자체가 없다(diploid 단일 index 라 리드는 한 번만 정렬된다).
  따라서 기본값을 바꾸는 근거를 AlleleSeq2 에서 가져올 수 없다 — 별도 판단으로 남긴다.

★ 2026-09-17 수정 — AlleleSeq2 대조 반영
  ① duplicate 제외. AlleleSeq2 는 `RM_DUPLICATE_READS := on` 이 기본이고
     `FINAL_ALIGNMENT_FILENAME` 이 `rmdup.bam` 을 가리켜 **카운트 전에 물리적으로 제거**한다.
     우리는 17 이 `REMOVE_DUPLICATES=false` 로 flag 만 남기므로 여기서 거른다.
     dense coverage(18)는 GTEx 기준으로 duplicate 를 포함한다 — coverage 와 ASE count 의 duplicate 정책은
     의도적으로 다르다(coverage = GTEx, ASE = Castel et al. 2015 의 ASE QC).
  ② multimapper 보정 산출물 추가(`.mmap_adjusted.tsv`).
     AlleleSeq2 는 uniq/mmap 정렬을 별도 BAM·mpileup 으로 뽑아
     `filter_sites_w_mmaps.py adjust` 로 **약한 allele 쪽에 mmap 카운트를 더해** 불균형이
     multimapping 탓인지 본다. MAPQ 로 버리기만 하면 카운트는 보수적이 되지만 편향은 안 막힌다 —
     강한 allele 쪽 자리만 uniquely mappable 하면 불균형이 그대로 남는다.
     여기서는 gene 단위로 같은 구조를 만든다: uniq(MAPQ>=min) 과 mmap(0<MAPQ<min) 을
     서로소 집합으로 각각 배정한 뒤, 약한 쪽에 mmap 을 더해 방향이 뒤집히는지 표시한다.

★ hapA=GT[0] 규약
  HAP1 = `bcftools consensus -H 1` = GT[0] 이므로 aCount=hap1 로 두면 기존 ASE 진리값과 규약이
  일치한다(ase_validation_260820.md §3 의 방향 검증과 같은 결론).
  ⚠ 단 GT 의 `1|0`/`0|1` 방향은 **phase set 마다 독립**이라, HAP1 은 염색체 haplotype 이 아니라
  phase block 모자이크다. 패널 기반 진리값과 **방향 비교는 불가**하고 |ratio−0.5| 크기만 비교된다.
  근거: 같은 문서 §4. 유전자의 99.2% 는 단일 block 안에 있으므로 유전자 단위 측정 자체는 깨끗하다.
"""
from __future__ import annotations
import argparse, collections, hashlib, json, math, os, sys

import pysam
from rna_allelic_common import (
    MATES, EXON_BIN, split_bucket, load_exons, gene_set, restrict_to_common,
    genes_for, build_exon_index, genes_for_indexed, genes_of, assign_fragments,
)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rna_profile  # noqa: E402

def read_as(bam_path: str, mapq_min: int, include_duplicates: bool = False):
    """primary 정렬을 MAPQ 로 갈라 두 벌로 읽는다.

    반환 (uniq, mmap). 각각 query_name → [mate0, mate1], 원소는 (AS, contig, blocks) 또는 None.
      uniq: MAPQ >= mapq_min        — 기존 카운트 경로
      mmap: 0 < MAPQ <  mapq_min    — multimapper 보정용(AlleleSeq2 의 mmapreads BAM 대응)
    두 집합은 서로소다. MAPQ 0(정렬 위치 불확정)은 양쪽 모두에서 제외한다.

    duplicate 는 기본적으로 제외한다 — AlleleSeq2 는 카운트 전에 rmdup 으로 제거한다.
    어느 fragment 가 duplicate 인지는 호출부가 두 HAP 의 합집합으로 다시 정한다
    (drop_duplicate_fragments). 여기서는 이 BAM 에서 DUP flag 가 붙은 이름만 모아 돌려준다.
    """
    uniq: dict = {}
    mmap: dict = {}
    n_dup_skipped = 0
    dup_names: set = set()
    with pysam.AlignmentFile(bam_path, "rb") as bf:
        for r in bf.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if r.is_duplicate:
                n_dup_skipped += 1
                dup_names.add(r.query_name)
                if not include_duplicates:
                    continue
            mq = r.mapping_quality
            if mq >= mapq_min:
                target = uniq
            elif mq > 0:
                target = mmap
            else:
                continue
            try:
                a = int(r.get_tag("AS"))
            except KeyError:
                continue
            i = 0 if r.is_read1 else 1
            slot = target.get(r.query_name)
            if slot is None:
                slot = target[r.query_name] = [None, None]
            prev = slot[i]
            if prev is None or a > prev[0]:
                slot[i] = (a, r.reference_name, r.get_blocks())
    return uniq, mmap, n_dup_skipped, dup_names


STAR_GENOMIC_LENGTH_LOG2SCALE = rna_profile.STAR_GENOMIC_LENGTH_LOG2SCALE   # 16 이 쓰는 값(STAR·GTEx 기본)


def star_length_term(span: int, scale: float) -> int:
    """STAR 가 AS 에 더한 genomic-length 항 (source/stitchWindowAligns.cpp 221–225):
    int(ceil(log2(span) * scoreGenomicLengthLog2scale - 0.5)). scale == 0 이면 STAR 도 더하지 않는다."""
    if scale == 0 or span <= 0:
        return 0
    return int(math.ceil(math.log2(span) * scale - 0.5))


def remove_length_term(table: dict, scale: float) -> int:
    """HAP 비교 전에 AS 에서 STAR 의 genomic-length 항을 뺀다 (16 은 GTEx 기본 -0.25 로 정렬).

    ★ 왜
      STAR 의 AS 는 쌍 점수이고 여기에 ceil(log2(span) × s − 0.5) 가 들어 있다(span = 두 mate 를 합친 genomic 길이).
      s = −0.25 이면 span 이 4·64·1,024·16,384·262,144 bp 를 넘을 때마다 1 점씩 깎인다. 같은 fragment 의 HAP1·HAP2
      span 은 그 안의 indel 만큼 다르므로, 두 span 이 경계를 사이에 두면 서열 일치와 무관하게 span 이 짧은 HAP 이
      1 점 이긴다. main 정렬(coverage)은 GTEx 그대로 두고, 여기서 이 항만 정확히 되돌려 가짜 HAP 차이를 없앤다.
      같은 contig 의 mate 끼리만 span 을 잰다(다르면 STAR 도 한 정렬로 잇지 않는다 → 보정하지 않음).
    반환: 보정한 fragment 수.
    """
    if scale == 0:
        return 0
    n = 0
    for slot in table.values():
        mates = [m for m in slot if m is not None]
        if not mates or len({m[1] for m in mates}) != 1:
            continue
        lo = min(b[0][0] for m in mates for b in [m[2]] if b)
        hi = max(b[-1][1] for m in mates for b in [m[2]] if b)
        term = star_length_term(hi - lo, scale)
        if term:
            for i, m in enumerate(slot):
                if m is not None:
                    slot[i] = (m[0] - term, m[1], m[2])
            n += 1
    return n


def drop_duplicate_fragments(dup_names: set, *tables: dict) -> int:
    """두 HAP 중 어느 한쪽에서라도 duplicate 인 fragment 를 모든 표에서 뺀다. 뺀 수를 돌려준다.

    ★ 왜 (2026-09-28)
      14 는 MarkDuplicates 를 HAP1·HAP2 BAM 에 따로 돌린다. duplicate 판정은 정렬 위치
      (unclipped 5' 끝)로 묶이므로, HAP 사이 indel·soft-clip 차이로 같은 fragment 가 한쪽에서만
      duplicate 가 될 수 있다. 그러면 그 fragment 는 한쪽에서만 빠져 only1/only2 로 배정되고,
      기술적 인공물이 allele 증거가 된다. NA20503 파일럿 실측: 두 HAP 에 모두 정렬된
      fragment 261만 개 중 11,790개(HAP1 에서만 9,589 / HAP2 에서만 2,201)가 이런 경우였다.
      Castel et al. 2015(Genome Biol 16:195)는 AE 분석에서 duplicate 제거를 기본으로 권하되
      제거가 allele 에 편향되지 않아야 한다고 한다. fragment 단위로 합집합을 취하면 두 HAP 에
      같은 fragment 집합이 남는다.
    """
    removed = 0
    for table in tables:
        for name in dup_names:
            if table.pop(name, None) is not None:
                removed += 1
    return removed


# exon 탐색과 구번호 fragment 배정은 rna_allelic_common의 호환 구현을 사용한다.
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam1", required=True); ap.add_argument("--bam2", required=True)
    ap.add_argument("--exons1", required=True); ap.add_argument("--exons2", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-mapq", type=int, default=255, help="STAR 는 unique 정렬에 255 를 준다")
    ap.add_argument("--require-both-haps", action="store_true",
                    help="한쪽 haplotype 에만 정렬된 fragment 를 제외한다(민감도 분석용). "
                         "기본값은 제외하지 않음 — 기존 동작 유지")
    ap.add_argument("--include-duplicates", action="store_true",
                    help="Picard DUP flag 가 붙은 정렬을 카운트에 포함한다. 기본값은 제외이며 "
                         "ASE QC 목적(Castel et al. 2015). 18 coverage 는 GTEx 기준으로 duplicate 를 포함한다")
    a = ap.parse_args()

    A, A_mm, dup1, dupnames1 = read_as(a.bam1, a.min_mapq, a.include_duplicates)
    B, B_mm, dup2, dupnames2 = read_as(a.bam2, a.min_mapq, a.include_duplicates)
    dup_union = dupnames1 | dupnames2
    dup_inconsistent = len(dupnames1 ^ dupnames2)
    n_dup_union_removed = 0
    if not a.include_duplicates:
        n_dup_union_removed = drop_duplicate_fragments(dup_union, A, A_mm, B, B_mm)
    n_len_corrected = sum(remove_length_term(t, STAR_GENOMIC_LENGTH_LOG2SCALE) for t in (A, A_mm, B, B_mm))
    raw1, raw2 = load_exons(a.exons1), load_exons(a.exons2)
    raw1, raw2, n_gene_dropped = restrict_to_common(raw1, raw2)
    ex1 = build_exon_index(raw1)
    ex2 = build_exon_index(raw2)
    if n_gene_dropped:
        print(f"[allelic] 한쪽 hap 에만 있는 유전자 {n_gene_dropped:,}개 제외 "
              f"(비교 불가 — 구조적 편향 방지)")

    cnt, st = assign_fragments(A, B, ex1, ex2, a.require_both_haps)
    aC, bC = cnt["a"], cnt["b"]
    both_a, both_b = cnt["both_a"], cnt["both_b"]
    split_a, split_b = cnt["split_a"], cnt["split_b"]
    n_tie, n_only1, n_only2 = st["tie"], st["only1"], st["only2"]
    n_assigned, n_nocommon, n_dropped_single = st["assigned"], st["nocommon"], st["dropped_single"]

    # multimapper 보정 — 같은 규칙을 mmap 집합에 적용한다(서로소이므로 이중 계수 없음).
    mm_cnt, mm_st = assign_fragments(A_mm, B_mm, ex1, ex2, a.require_both_haps)
    mmA, mmB = mm_cnt["a"], mm_cnt["b"]

    with open(a.out, "w") as fh:
        fh.write("gene_id\taCount\tbCount\ttotalCount\tratio\n")
        for g in sorted(set(aC) | set(bC)):
            x, y = aC[g], bC[g]; t = x + y
            fh.write(f"{g}\t{x}\t{y}\t{t}\t{x/t if t else 0.5:.6f}\n")
    split_path = os.path.splitext(a.out)[0] + ".split_half.tsv"
    with open(split_path, "w") as fh:
        fh.write("gene_id\tsplit0_aCount\tsplit0_bCount\tsplit0_totalCount\tsplit0_ratio\t"
                 "split1_aCount\tsplit1_bCount\tsplit1_totalCount\tsplit1_ratio\n")
        for g in sorted(set(aC) | set(bC)):
            a0, b0 = split_a[0][g], split_b[0][g]
            a1, b1 = split_a[1][g], split_b[1][g]
            t0, t1 = a0 + b0, a1 + b1
            r0 = a0 / t0 if t0 else 0.5
            r1 = a1 / t1 if t1 else 0.5
            fh.write(f"{g}\t{a0}\t{b0}\t{t0}\t{r0:.6f}\t"
                     f"{a1}\t{b1}\t{t1}\t{r1:.6f}\n")
    both_path = os.path.splitext(a.out)[0] + ".both_haps.tsv"
    with open(both_path, "w") as fh:
        fh.write("gene_id\taCount\tbCount\ttotalCount\tratio\n")
        for g in sorted(set(both_a) | set(both_b)):
            x, y = both_a[g], both_b[g]; t = x + y
            fh.write(f"{g}\t{x}\t{y}\t{t}\t{x/t if t else 0.5:.6f}\n")

    # ── multimapper 보정 (AlleleSeq2 filter_sites_w_mmaps.py adjust 대응) ──
    #   약한 쪽에 mmap 카운트를 몰아준 뒤에도 방향이 유지되는지 본다. 뒤집히면
    #   그 유전자의 불균형은 uniquely mappable 한 자리 편중으로 설명될 수 있다는 뜻이다.
    mmap_path = os.path.splitext(a.out)[0] + ".mmap_adjusted.tsv"
    with open(mmap_path, "w") as fh:
        fh.write("gene_id\taCount\tbCount\taMmap\tbMmap\t"
                 "adj_aCount\tadj_bCount\tadj_ratio\tmmap_survives\n")
        for g in sorted(set(aC) | set(bC)):
            x, y = aC[g], bC[g]
            mx, my = mmA[g], mmB[g]
            if x == y:
                ax, ay = x + mx, y + my          # 방향이 없으면 그대로 더한다
            elif x > y:
                ax, ay = x, y + mx + my          # 약한 쪽(b)에 mmap 전량을 준다
            else:
                ax, ay = x + mx + my, y
            at = ax + ay
            # 원래 방향이 보정 후에도 같은 쪽으로 남는가
            survives = (x > y and ax > ay) or (x < y and ax < ay)
            fh.write(f"{g}\t{x}\t{y}\t{mx}\t{my}\t{ax}\t{ay}\t"
                     f"{ax/at if at else 0.5:.6f}\t{1 if survives else 0}\n")
    tot = len(set(A) | set(B))

    # sidecar — TSV 컬럼은 그대로 두어야 하위 소비자가 안 깨지므로, 단위와 조건은 옆에 적는다.
    # 이게 없으면 "이 파일이 read 단위인가 fragment 단위인가"를 파일만 보고 알 수 없다.
    meta_path = os.path.splitext(a.out)[0] + ".meta.json"
    with open(meta_path, "w") as fh:
        json.dump({
            "schema": "mfm_allelic_counts_meta_v4",
            "counting_unit": "fragment",
            "script_revision": "2026-09-28_hap_union_duplicates",
            "duplicate_unit": "fragment, union of both HAP DUP flags",
            "duplicate_fragments_union": len(dup_union),
            "duplicate_fragments_flag_inconsistent_between_haps": dup_inconsistent,
            "duplicate_union_entries_removed": n_dup_union_removed,
            "min_mapq": a.min_mapq,
            "as_comparison": f"STAR AS minus genomic-length term ceil(log2(span)*s-0.5), s={STAR_GENOMIC_LENGTH_LOG2SCALE:g}",
            "fragments_as_length_term_removed": n_len_corrected,
            "require_both_haps": bool(a.require_both_haps),
            "include_duplicates": bool(a.include_duplicates),
            "duplicate_policy": ("included (--include-duplicates)" if a.include_duplicates
                                 else "excluded (BAM DUP flag, 두 HAP 합집합; ASE QC 목적. 18 coverage 는 duplicate 포함)"),
            "duplicate_primary_records_seen": {"bam1": dup1, "bam2": dup2},
            "bam1": a.bam1, "bam2": a.bam2,
            "fragments_total": tot,
            "fragments_assigned": n_assigned,
            "fragments_tie": n_tie,
            "fragments_no_common_mate": n_nocommon,
            "fragments_only_hap1": n_only1,
            "fragments_only_hap2": n_only2,
            "fragments_dropped_single_hap": n_dropped_single,
            "genes": len(set(aC) | set(bC)),
            "genes_excluded_hap_asymmetric": n_gene_dropped,
            "split_half_path": split_path,
            "split_half_method": "blake2b(query_name, person=mfm-ase), digest_byte0 & 1",
            "both_haps_only_path": both_path,
            "mmap_adjusted_path": mmap_path,
            "mmap_adjustment": {
                "method": "AlleleSeq2 filter_sites_w_mmaps.py adjust의 gene 단위 대응",
                "uniq_mapq": f">= {a.min_mapq}",
                "mmap_mapq": f"0 < MAPQ < {a.min_mapq}",
                "fragments_assigned": mm_st["assigned"],
                "fragments_tie": mm_st["tie"],
                "genes": len(set(mmA) | set(mmB)),
            },
            "note": ("2026-09-04 이전 산출물은 read 단위라 카운트가 약 2배이고 이항검정이 "
                     "과도하게 유의하다. 이 sidecar 가 없는 allelic_counts.tsv 는 read 단위로 간주할 것."),
        }, fh, ensure_ascii=False, indent=2)
        fh.write("\n")
    print(f"[allelic] fragment {tot:,} / 배정 {n_assigned:,} ({100*n_assigned/max(tot,1):.1f}%) / "
          f"동점 제외 {n_tie:,} / 공통 mate 없음 {n_nocommon:,} / "
          f"한쪽만 정렬 H1 {n_only1:,} H2 {n_only2:,}"
          + (f" / --require-both-haps 로 제외 {n_dropped_single:,}" if a.require_both_haps else ""))
    print(f"[allelic] 단위 = fragment (2026-09-04 이전 산출물은 read 단위라 카운트가 약 2배다)")
    print(f"[allelic] duplicate {'포함' if a.include_duplicates else '제외'} "
          f"(primary DUP 레코드 H1 {dup1:,} / H2 {dup2:,}; fragment 합집합 {len(dup_union):,}, "
          f"HAP 간 불일치 {dup_inconsistent:,})")
    print(f"[allelic] mmap 보정: fragment {mm_st['assigned']:,} 배정 → {mmap_path}")
    print(f"[allelic] 유전자 {len(set(aC)|set(bC)):,} → {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
