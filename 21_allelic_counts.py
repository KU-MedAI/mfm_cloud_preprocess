#!/usr/bin/env python3
"""두 haplotype 정렬을 **fragment 별 AS 비교**로 합쳐 대립 카운트를 만든다 (ASE 진리값).

    python 21_allelic_counts.py --bam1 <HAP1.bam> --bam2 <HAP2.bam> \
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

★ 2026-10-02 수정 — 검증된 HAP variant 를 덮는 fragment 만 allele 증거로 인정 (20 단계 산출물 사용)
  위 설계는 "두 HAP 정렬이 다르면 그 fragment 가 HAP 사이 차이를 지났다"는 전제에 기댄다. HG00096 EDA 에서
  이 전제가 깨지는 두 경우가 나왔다.
    ① 한쪽 HAP 에서만 unique(only1/only2): 다른 곳 pseudogene 의 변이로 유일성이 바뀐다. RPS15 3:707 은 gene 위에
       HAP 차이가 0 개였고 같은 read 의 AS 가 두 HAP 에서 같았다. 이진 ASE 적격 211 중 82, 유의 379 중 218 gene 이
       신호의 절반 이상을 이 경로에서 얻었다.
    ② VCF 오류: 실제로 없는 variant 가 HAP 서열에 들어가 점수를 깎는다(PKN1 chr19:14,469,962 TC>T, HiFi 0/35).
  그래서 배정 전에 fragment 정렬 블록(해당 HAP 좌표)이 20 의 hap_variants.tsv 구간과 겹치는지 본다.
    - validated 가 아닌 variant(unsupported·not_applied·low_coverage)를 어느 HAP 정렬에서든 덮으면 → 버림
    - validated variant 를 이긴 HAP(only 경로는 그 HAP) 정렬에서 덮지 않으면 → 버림(HAP 차이 밖에서 생긴 신호)
    - 둘 다 아니면 기존 규칙대로 배정
  fragment 한 번만 세기·AS 비교·길이 항 제거·duplicate 합집합 제외는 그대로다.

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
    MATES, split_bucket, load_exons, gene_set, genes_for, build_exon_index,
    genes_for_indexed, genes_of,
)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rna_profile  # noqa: E402

def sha256_file(path: str, chunk: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


class HapEvidence:
    """20 단계 hap_variants.tsv → HAP 좌표 구간 색인. validated 구간과 그 외(근거로 못 쓰는) 구간을 따로 둔다."""

    def __init__(self, path: str):
        import bisect  # noqa: F401
        good: dict = collections.defaultdict(list); bad: dict = collections.defaultdict(list)
        self.counts = collections.Counter()
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            col = {k: i for i, k in enumerate(head)}
            for k in ("contig", "status", "kind", "h1_start", "h1_end", "h2_start", "h2_end"):
                if k not in col:
                    raise SystemExit(f"[allelic] hap_variants 열 없음: {k} ({path})")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                st = f[col["status"]]; self.counts[st] += 1
                pad = 1 if f[col["kind"]] == "INDEL" else 0      # 결실 접합부를 지나는 정렬도 겹침으로 본다
                for h in (1, 2):
                    iv = (int(f[col[f"h{h}_start"]]), int(f[col[f"h{h}_end"]]) + pad)
                    (good if st == "validated" else bad)[(h, f[col["contig"]])].append(iv)
        self.idx = {}
        for name, d in (("good", good), ("bad", bad)):
            for key, ivs in d.items():
                ivs.sort()
                starts = [a for a, _ in ivs]; ends = [b for _, b in ivs]
                maxlen = max(b - a for a, b in ivs)
                self.idx[(name, *key)] = (starts, ends, maxlen)

    def _hit(self, name: str, h: int, slot) -> bool:
        import bisect
        for m in slot:
            if m is None:
                continue
            e = self.idx.get((name, h, m[1]))
            if e is None:
                continue
            starts, ends, maxlen = e
            for bs, be in m[2]:
                i = bisect.bisect_left(starts, bs - maxlen); j = bisect.bisect_left(starts, be)
                if any(ends[k] > bs for k in range(i, j)):
                    return True
        return False

    def validated(self, h: int, slot) -> bool:
        return self._hit("good", h, slot)

    def unvalidated(self, h: int, slot) -> bool:
        return self._hit("bad", h, slot)


# 순차 API를 보존하며 병렬 입력 엔진은 공통 모듈에서 제공함.
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parent/'src'))
from integrate_wgs_rna.allelic_input import read_as, read_pair_as
from integrate_wgs_rna.allelic_parallel import assign_partitioned


STAR_GENOMIC_LENGTH_LOG2SCALE = rna_profile.STAR_GENOMIC_LENGTH_LOG2SCALE   # 16 이 쓰는 값(STAR·GTEx 기본)


def star_length_term(span: int, scale: float) -> int:
    """STAR 가 AS 에 더한 genomic-length 항 (source/stitchWindowAligns.cpp 221–225):
    int(ceil(log2(span) * scoreGenomicLengthLog2scale - 0.5)). scale == 0 이면 STAR 도 더하지 않는다."""
    if scale == 0 or span <= 0:
        return 0
    return int(math.ceil(math.log2(span) * scale - 0.5))


def remove_length_term(table: dict, scale: float) -> int:
    """기존 STAR 길이 보정식을 두 mate의 직접 비교로 계산함."""
    if scale == 0:
        return 0
    corrected = 0
    for slot in table.values():
        left, right = slot
        if left is None and right is None:
            continue
        if left is not None and right is not None and left[1] != right[1]:
            continue
        a = left[2] if left is not None else ()
        b = right[2] if right is not None else ()
        if a and b:
            start = min(a[0][0], b[0][0])
            end = max(a[-1][1], b[-1][1])
        elif a:
            start, end = a[0][0], a[-1][1]
        elif b:
            start, end = b[0][0], b[-1][1]
        else:
            raise ValueError('min() arg is an empty sequence')
        term = star_length_term(end - start, scale)
        if term:
            if left is not None:
                slot[0] = (left[0] - term, left[1], left[2])
            if right is not None:
                slot[1] = (right[0] - term, right[1], right[2])
            corrected += 1
    return corrected


def drop_duplicate_fragments(dup_names: set, *tables: dict) -> int:
    """두 HAP 중 어느 한쪽에서라도 duplicate 인 fragment 를 모든 표에서 뺀다. 뺀 수를 돌려준다.

    ★ 왜 (2026-09-28)
      17은 MarkDuplicates를 HAP1·HAP2 BAM에 따로 돌린다. duplicate 판정은 정렬 위치
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


def restrict_to_common(ex1, ex2):
    """양쪽 hap 에 모두 존재하는 유전자만 남긴다. (ex1, ex2, 버린 유전자 수) 반환.

    ★ 왜 필요한가 (2026-09-17)
      배정된 fragment 는 **이긴 haplotype 의 exon 목록**으로만 유전자를 찾는다
      (HAP1 승리 → ex1, HAP2 승리 → ex2). 그래서 어떤 유전자가 한쪽 exon 목록에만 있으면
      그쪽으로만 카운트가 쌓이고 반대쪽은 **구조적으로 0** 이 된다.
      → 그 유전자는 "완벽한 편향 ASE" 로 보이지만 전부 인공물이다.

      Stage 14의 기본값이 snap이던 시절에는 lift가 아무것도 버리지 않아 두 목록이 항상 같았고
      (실측: ONT 10샘플 전부 55,746 = 55,746) 이 문제가 드러나지 않았다.
      drop 으로 바꾸면 lift 실패 위치가 hap 마다 달라 반드시 어긋난다.
      한쪽에서 위치를 못 찾는 유전자는 haplotype 간 비교 자체가 불가능하므로 교집합만 센다.
    """
    common = gene_set(ex1) & gene_set(ex2)
    dropped = len(gene_set(ex1) | gene_set(ex2)) - len(common)
    out = []
    for ex in (ex1, ex2):
        filt = collections.defaultdict(list)
        for c, items in ex.items():
            keep = [(s, e, g) for s, e, g in items if g in common]
            if keep:
                filt[c] = keep
        out.append(filt)
    return out[0], out[1], dropped


def assign_fragments(A, B, ex1, ex2, require_both_haps: bool, ev: "HapEvidence", *, keys=None):
    """fragment 별 AS 비교 → gene 카운트. 배정 규칙 전체가 여기 모여 있다.

    uniq 집합과 mmap 집합에 같은 규칙을 적용하기 위해 함수로 분리했다.
    """
    aC = collections.Counter(); bC = collections.Counter()
    both_a = collections.Counter(); both_b = collections.Counter()
    split_a = (collections.Counter(), collections.Counter())
    split_b = (collections.Counter(), collections.Counter())
    st = dict(tie=0, only1=0, only2=0, assigned=0, nocommon=0, dropped_single=0,
              drop_unvalidated_variant=0, drop_no_validated_variant=0)
    for key in (set(A) | set(B)) if keys is None else keys:
        ra, rb = A.get(key), B.get(key)
        if ra is not None and rb is not None:
            # 양쪽에 다 정렬된 mate 에서만 점수를 합산해 비교한다.
            common = [i for i in MATES if ra[i] is not None and rb[i] is not None]
            if not common:
                # 양쪽에 있긴 한데 겹치는 mate 가 없다(mate1 은 H1 에만, mate2 는 H2 에만).
                # 공정한 비교가 불가능하므로 무정보로 버린다.
                st["nocommon"] += 1; continue
            s1 = sum(ra[i][0] for i in common)
            s2 = sum(rb[i][0] for i in common)
            if s1 == s2:
                st["tie"] += 1; continue
            win, slot = (1, ra) if s1 > s2 else (2, rb)
            if ev.unvalidated(1, ra) or ev.unvalidated(2, rb):
                st["drop_unvalidated_variant"] += 1; continue
            if not ev.validated(win, slot):
                st["drop_no_validated_variant"] += 1; continue
            bucket = split_bucket(key)
            if s1 > s2:
                for g in genes_of(ex1, ra):
                    aC[g] += 1; split_a[bucket][g] += 1; both_a[g] += 1
            else:
                for g in genes_of(ex2, rb):
                    bC[g] += 1; split_b[bucket][g] += 1; both_b[g] += 1
        else:
            if require_both_haps:
                st["dropped_single"] += 1; continue
            h, slot = (1, ra) if ra is not None else (2, rb)
            if ev.unvalidated(h, slot):
                st["drop_unvalidated_variant"] += 1; continue
            if not ev.validated(h, slot):
                st["drop_no_validated_variant"] += 1; continue
            bucket = split_bucket(key)
            if ra is not None:
                st["only1"] += 1
                for g in genes_of(ex1, ra):
                    aC[g] += 1; split_a[bucket][g] += 1
            else:
                st["only2"] += 1
                for g in genes_of(ex2, rb):
                    bC[g] += 1; split_b[bucket][g] += 1
        st["assigned"] += 1
    return dict(a=aC, b=bC, both_a=both_a, both_b=both_b,
                split_a=split_a, split_b=split_b), st


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam1", required=True); ap.add_argument("--bam2", required=True)
    ap.add_argument("--exons1", required=True); ap.add_argument("--exons2", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--hap-variants", required=True,
                    help="20_hap_variants.py 의 hap_variants.tsv — validated variant 를 덮는 fragment 만 allele 증거로 쓴다")
    ap.add_argument("--min-mapq", type=int, default=255, help="STAR 는 unique 정렬에 255 를 준다")
    ap.add_argument("--require-both-haps", action="store_true",
                    help="한쪽 haplotype 에만 정렬된 fragment 를 제외한다(민감도 분석용). "
                         "기본값은 제외하지 않음 — 기존 동작 유지")
    ap.add_argument("--include-duplicates", action="store_true",
                    help="Picard DUP flag 가 붙은 정렬을 카운트에 포함한다. 기본값은 제외이며 "
                         "ASE QC 목적(Castel et al. 2015). 18 coverage 는 GTEx 기준으로 duplicate 를 포함한다")
    a = ap.parse_args()

    jobs=os.environ.get('RNA_ALLELIC_JOBS','auto')
    workers=int(os.environ.get('MFM_SAMPLE_CPU_BUDGET','1')) if jobs=='auto' else int(jobs)
    if jobs=='auto' and sum(os.path.getsize(p) for p in [a.bam1,a.bam2])<64*1024**2:workers=1
    raw1, raw2 = load_exons(a.exons1), load_exons(a.exons2)
    raw1, raw2, n_gene_dropped = restrict_to_common(raw1, raw2)
    ex1 = build_exon_index(raw1)
    ex2 = build_exon_index(raw2)
    if n_gene_dropped:
        print(f"[allelic] 한쪽 hap 에만 있는 유전자 {n_gene_dropped:,}개 제외 "
              f"(비교 불가 — 구조적 편향 방지)")

    ev = HapEvidence(a.hap_variants)
    from integrate_wgs_rna.rna_ase_shards import selected,assign_from_bams
    from integrate_wgs_rna.gpu_compare_resources import cpu_limit
    workers=min(max(1,workers),cpu_limit(int(os.environ.get('MFM_SAMPLE_CPU_BUDGET',str(workers)))))
    from integrate_wgs_rna.rna_wgs_native import cpu_upper
    upper=cpu_upper(workers)
    if jobs!='auto':upper=min(upper,int(jobs))
    partitioned=selected([a.bam1,a.bam2],upper)
    if partitioned:
        cnt,st,mm_cnt,mm_st,au,input_execution,assignment_execution=assign_from_bams(
            [a.bam1,a.bam2],assign_fragments,drop_duplicate_fragments,remove_length_term,
            ex1,ex2,a.require_both_haps,ev,a.min_mapq,a.include_duplicates,STAR_GENOMIC_LENGTH_LOG2SCALE,upper,initial_workers=workers,native_policy='validated_hap_v1')
        if os.environ.get('RNA_ALLELIC_BENCHMARK','0')=='1':
            from integrate_wgs_rna.rna_ase_shards import verify_legacy
            check=verify_legacy([a.bam1,a.bam2],assign_fragments,drop_duplicate_fragments,remove_length_term,
                ex1,ex2,a.require_both_haps,ev,a.min_mapq,a.include_duplicates,STAR_GENOMIC_LENGTH_LOG2SCALE,(cnt,st,mm_cnt,mm_st,au))
            input_execution.update(serial_equal=True,legacy_seconds=check['input_seconds'],serial_equal_scope=check['scope'])
            input_execution['input_read_speedup']=check['input_seconds']/max(.000001,input_execution['seconds'])
            assignment_execution.update(serial_equal=True,serial_seconds=check['assignment_seconds'],serial_equal_scope=check['scope'])
            assignment_execution['speedup']=check['assignment_seconds']/max(.000001,assignment_execution['seconds'])
        dup1,dup2,dup_union_n,dup_inconsistent=au['dup1'],au['dup2'],au['dup_union'],au['dup_inconsistent']
        n_dup_union_removed,n_len_corrected,tot=au['removed'],au['corrected'],au['total']
    else:
        pair,input_execution=read_pair_as(a.bam1,a.bam2,a.min_mapq,a.include_duplicates,
            workers=max(1,workers),benchmark=os.environ.get('RNA_ALLELIC_BENCHMARK','0')=='1',
            scratch=str(__import__('pathlib').Path(a.out).resolve().parent))
        (A,A_mm,dup1,dupnames1),(B,B_mm,dup2,dupnames2)=pair
        dup_union = dupnames1 | dupnames2
        dup_inconsistent = len(dupnames1 ^ dupnames2)
        n_dup_union_removed = 0
        if not a.include_duplicates:
            n_dup_union_removed = drop_duplicate_fragments(dup_union, A, A_mm, B, B_mm)
        n_len_corrected = sum(remove_length_term(t, STAR_GENOMIC_LENGTH_LOG2SCALE) for t in (A, A_mm, B, B_mm))
        (cnt,st),assignment_execution=assign_partitioned(assign_fragments,A,B,ex1,ex2,a.require_both_haps,ev,
            workers=max(1,workers),elastic=jobs=='auto',benchmark=os.environ.get('RNA_ALLELIC_BENCHMARK','0')=='1',
            batch_size=int(os.environ.get('RNA_ALLELIC_BATCH_SIZE','50000')))
        mm_cnt,mm_st=assign_fragments(A_mm,B_mm,ex1,ex2,a.require_both_haps,ev)
        dup_union_n=len(dup_union);tot=len(set(A)|set(B))
    aC, bC = cnt["a"], cnt["b"]
    both_a, both_b = cnt["both_a"], cnt["both_b"]
    split_a, split_b = cnt["split_a"], cnt["split_b"]
    n_tie, n_only1, n_only2 = st["tie"], st["only1"], st["only2"]
    n_assigned, n_nocommon, n_dropped_single = st["assigned"], st["nocommon"], st["dropped_single"]

    # multimapper 보정 — 같은 규칙을 mmap 집합에 적용한다(서로소이므로 이중 계수 없음).
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

    # sidecar — TSV 컬럼은 그대로 두어야 하위 소비자가 안 깨지므로, 단위와 조건은 옆에 적는다.
    # 이게 없으면 "이 파일이 read 단위인가 fragment 단위인가"를 파일만 보고 알 수 없다.
    meta_path = os.path.splitext(a.out)[0] + ".meta.json"
    with open(meta_path, "w") as fh:
        json.dump({
            "schema": "mfm_allelic_counts_meta_v5",
            "resource_execution": {"input":input_execution,"assignment":assignment_execution,
                "serial_equal":(bool(input_execution.get("serial_equal")) and bool(assignment_execution.get("serial_equal"))) if os.environ.get("RNA_ALLELIC_BENCHMARK","0")=="1" else None,
                "assignment_speedup":assignment_execution.get("speedup") if assignment_execution.get("parallel_observed") else 1.0,
                "assignment_parallel_observed":bool(assignment_execution.get("parallel_observed")),"input_read_speedup":input_execution.get("input_read_speedup")},
            "counting_unit": "fragment",
            "script_revision": "2026-10-02_validated_hap_variant_evidence",
            "evidence_rule": ("fragment 가 20 의 validated HAP variant 를 이긴 HAP(only 경로는 그 HAP) 정렬에서 덮을 때만 배정. "
                              "validated 가 아닌 variant 를 어느 HAP 정렬에서든 덮으면 버림"),
            "hap_variants": a.hap_variants,
            "hap_variants_sha256": sha256_file(a.hap_variants),
            "hap_variants_status_counts": dict(ev.counts),
            "fragments_dropped_unvalidated_variant": st["drop_unvalidated_variant"],
            "fragments_dropped_no_validated_variant": st["drop_no_validated_variant"],
            "duplicate_unit": "fragment, union of both HAP DUP flags",
            "duplicate_fragments_union": dup_union_n,
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
            # Stage 23/final validator가 서로 다른 실행의 TSV를 조용히 섞지 않게
            # 실제 소비 파일 내용을 묶는다. 기존 TSV schema는 바꾸지 않는다.
            "output_sha256": {
                "counts": sha256_file(a.out),
                "split_half": sha256_file(split_path),
                "both_haps": sha256_file(both_path),
                "mmap_adjusted": sha256_file(mmap_path),
            },
            "mmap_adjustment": {
                "method": "AlleleSeq2 filter_sites_w_mmaps.py adjust의 gene 단위 대응",
                "uniq_mapq": f">= {a.min_mapq}",
                "mmap_mapq": f"0 < MAPQ < {a.min_mapq}",
                "fragments_assigned": mm_st["assigned"],
                "fragments_tie": mm_st["tie"],
                "fragments_dropped_unvalidated_variant": mm_st["drop_unvalidated_variant"],
                "fragments_dropped_no_validated_variant": mm_st["drop_no_validated_variant"],
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
    print(f"[allelic] 근거 규칙: 검증 안 된 variant 위 제외 {st['drop_unvalidated_variant']:,} / "
          f"검증된 HAP variant 밖 제외 {st['drop_no_validated_variant']:,}")
    print(f"[allelic] 단위 = fragment (2026-09-04 이전 산출물은 read 단위라 카운트가 약 2배다)")
    print(f"[allelic] duplicate {'포함' if a.include_duplicates else '제외'} "
          f"(primary DUP 레코드 H1 {dup1:,} / H2 {dup2:,}; fragment 합집합 {dup_union_n:,}, "
          f"HAP 간 불일치 {dup_inconsistent:,})")
    print(f"[allelic] mmap 보정: fragment {mm_st['assigned']:,} 배정 → {mmap_path}")
    print(f"[allelic] 유전자 {len(set(aC)|set(bC)):,} → {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
