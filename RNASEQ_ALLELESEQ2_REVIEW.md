> **이식 기록 (2026-09-25):** 이 문서는 post-training 참고 디렉토리에서 작성된 2026-09-17 검토 기록이다. 당시의 “원본 동결” 문구와 경로 설명은 역사적 맥락이며, 승인 후 RNA 13~19 구현을 현재 v2에 선택 이식했다.

# RNA-seq 파이프라인 코드 리뷰 — 변동 내역 및 AlleleSeq2 대조

작성일: 2026-09-17 (같은 날 §2.2·§2.3 정정 및 §3 반영 완료)
대상 코드: `/NHNHOME/korea/post-training/processing_code/` (원본 `/NHNHOME/korea/mfm_data/process_code/v2/` 2026-09-16 시점 복사본)
대조 대상: [gersteinlab/AlleleSeq2](https://github.com/gersteinlab/AlleleSeq2) `1f3d424` (2021-04-07)

> 원본 `mfm_data/process_code/v2/`는 동결한다. 수정은 이 디렉토리에서만 한다.

**저장소 동일성 확인.** [trgaleev/AlleleSeq2](https://github.com/trgaleev/AlleleSeq2) `cfe8acf`는
`gersteinlab/master`를 머지한 fork이고 `.git`을 뺀 **전 파일이 바이트 단위로 동일**하다.
어느 쪽을 대조해도 결과가 같다.

**⚠ 이 문서의 초판에는 오류가 하나 있었다.** §2.2의 `sjdbScore` 항목이 AlleleSeq2의 RNA 경로를
잘못 인용했다. §2.3-② 아래 「정정」과 §3의 조치 2를 볼 것. 소스를 다시 읽고 확인한 결과다.

---

## 0. 리뷰 방법 / 기준선

`process_code/`는 git 저장소가 아니다. "변동"의 기준선은 다음 두 스냅샷을 사용했다.

| 기준선 | 경로 | 시점 |
|---|---|---|
| 1차 | `mfm_data/process_code_260905.tar.gz` | 2026-09-05 |
| 2차 | `mfm_data/process_code/v2_old/` | 2026-08-25 |

아래 "변동 내역"은 **9/5 스냅샷 → 9/16 현재** 기준이다.

AlleleSeq2는 `PIPELINE.mk`, `makePersonalGenome.mk`, `pileup2counts.py`, `filter_cnv_sites.py`,
`filter_phase_warnings.py`, `filter_sites_w_mmaps.py`, `get_hetSNV_bed.py`를 직접 읽고 대조했다.

---

## 1. 변동 내역 (2026-09-05 → 2026-09-16)

### 1.1 신규 파일

| 파일 | 줄 수 | 내용 |
|---|---|---|
| `02_preflight.sh` | 320 | WGS stage 실행 전 sample별 입력 계약 검사 (platform / decode ref / read type / VCF source) |
| `01_prepare_analysis_bam.sh` | 45 | raw long-read uBAM/FASTQ → analysis reference 좌표 BAM 표준화 |
| **`18_ase_qc.py`** | **246** | **ASE label 깊이 · split-half 재현성 · single-HAP 민감도 QC (RNA 신규 단계)** |
| `src/integrate_wgs_rna/prepare_alignment.py` | 325 | Stage 04 정렬 준비 |
| `src/integrate_wgs_rna/input_manifest.py` | 278 | 입력 manifest 검증 |
| `src/integrate_wgs_rna/public_vcf.py` | 358 | 공개 VCF 라우팅 |
| `src/integrate_wgs_rna/cohort_priority.py` | 365 | 코호트 우선순위 |
| `tests/`, `share_package/`, `cohorts/pacbio_1kgp/` | — | 테스트 스위트 · 배포 패키지 (전체 신규) |

### 1.2 수정 규모

```
config_v2.sh                +334 -41    manifest·코호트·플랫폼 프리셋 도입
src/.../personal_realign.py +313  -8
14_realign_rna_star.sh      +214 -29    ← RNA 핵심
run_pipeline.sh              +66  -8
src/.../cli.py               +68  -3
runner_local.sh              +68  -8
06_haplotag.sh               +54  -8
16_allelic_counts.py         +48  -7    ← RNA 핵심
03_prepare_vcf.sh            +36  -3
13_lift_annotation.py        +33  -4
slurm_rna.sbatch              +9  -4
15_coverage_zarr.py           +9  -2
run_rna_sample.sh            +17  -2
19_stamp_rna_metadata.py     +12  -0
```

### 1.3 RNA 경로 변경 상세

**`14_realign_rna_star.sh`**
1. BAM→FASTQ 단계에 pair 계약 검증 추가. 기존 `-0 /dev/null -s /dev/null`로 버리던 unpaired/singleton을
   실제 파일로 받아 레코드 수를 세고, R1 ≠ R2이면 `die`. `fastq_qc.json` 생성.
2. **Picard MarkDuplicates 단계 신설** — `REMOVE_DUPLICATES=false` (flag만 표기, 제거 안 함).
   `RNA_MARKDUP_REMOVE_DUPLICATES=0`만 허용하도록 가드.
3. contig 소스를 `$AUTOSOMES` → `$RNA_CONTIGS` (= `chroms_space $CHROMS`)로 변경.
4. skip 판정에 `mfm_v2_rna_realign_v2` 스키마 확인 추가 (구버전 산출물 재사용 방지).
5. `topmed_parameter_policy` 블록으로 TOPMed 대비 동일/상이 항목을 `rna_realign.json`에 명시.

**`16_allelic_counts.py`**
- **split-half 재현성 산출물 추가** — `blake2b(query_name, person=b"mfm-ase")` 최하위 비트로
  fragment를 0/1로 분할 → `.split_half.tsv`. 재실행·프로세스 무관하게 재현된다.
- **`.both_haps.tsv` 추가** — 양쪽 hap에 모두 정렬된 fragment만 집계한 민감도 분석용 산출물.
- meta schema `v1` → `v2`, `script_revision` = `2026-09-16_fragment_split_half`.

**`13_lift_annotation.py`**
- chain 파일 누락 시 조용히 건너뛰던 것을 `FileNotFoundError`로 변경.
- `--qc-json` 추가 (`mfm_gtf_lift_qc_v1`: 입력/출력 feature 수, snap 비율, exon 수).

**`15_coverage_zarr.py`**
- duplicate 정책 명시 — `samtools depth` 기본 필터로 DUP 제외, `--include-duplicates`로만 포함.

**`run_rna_sample.sh`**
- `18_ase_qc.py` 호출 추가, sidecar 3종(`split_half` / `both_haps` / `meta.json`) 존재 검증.

---

## 2. AlleleSeq2 대조

### 2.1 설계 차이 (의도된 분기)

| 항목 | AlleleSeq2 | v2 |
|---|---|---|
| 개인 유전체 | vcf2diploid → hap1+hap2를 **하나의 diploid FASTA로 합쳐 STAR index 1개** (+ unplaced/decoy/EBV contig) | bcftools consensus → **HAP1/HAP2 각각 별도 index 2개**, autosome만, decoy 없음 |
| 정렬 | 1회, `outFilterMultimapNmax 999999` → MAPQ 255 = **diploid 전체에서 unique** | 2회 독립, `outFilterMultimapNmax 20` → MAPQ 255 = **해당 haplotype 안에서만 unique** |
| 대립 배정 단위 | **hetSNV site** (lift된 좌표에서 mpileup base count) | **fragment의 AS(H1) vs AS(H2) 비교 → gene 단위** |
| 통계 | binomial + 시뮬레이션 FDR **AND** **beta-binomial (샘플별 overdispersion 추정, VGAM)** | exact binomial + BH FDR |
| 카운트 임계 | `Cntthresh_tot=6`, `Cntthresh_min=0` | `min_total_count=20` + split-half 안정성 |
| duplicate | Picard `REMOVE_DUPLICATES=true` (제거) | `false` (flag만) |
| 주석 lift | UCSC `liftOver -gff` + vcf2diploid chain, diploid GTF를 mapping 시점에 공급 | 자체 `13_lift_annotation.py` + `bcftools consensus --chain`, HAP별 GTF를 genomeGenerate 시점에 공급 |
| 영역 집계 | hetSNV 카운트를 bedtools로 영역에 집계 후 같은 통계 재적용 (`PIPELINE_aggregate_over_genomic_regions.mk`) | 처음부터 gene 단위 fragment 카운트 |

두 갈래 독립 정렬 후 AS를 비교하는 방식은 phASER/WASP 계열 접근이며 AlleleSeq2의 diploid 단일 index와
다른 정당한 설계다. 문제는 **AlleleSeq2가 diploid 단일 index를 쓰기 때문에 필요 없었던 보정들이
v2에는 빠져 있다**는 점이다.

### 2.2 STAR 파라미터 대조

| 옵션 | AlleleSeq2 (ASE) | v2 `31` | 비고 |
|---|---|---|---|
| `outFilterMismatchNoverReadLmax` | **0.03** | 0.1 (`NoverLmax`) | AlleleSeq2는 haplotype과 잘 안 맞는 리드를 아예 배제 |
| `outFilterMatchNminOverLread` | **0.95** | 0.33 | v2는 TOPMed 발현정량 기준 |
| `outFilterMultimapNmax` | 999999 | 20 | AlleleSeq2는 multimapper를 보존해 별도 보정 |
| `scoreGenomicLengthLog2scale` | **0.0** | 미지정 → STAR 기본 **−0.25** | §2.3 ② 참조 |
| `sjdbScore` | **ASB/ASCA(DNA)에서만 0**; ASE는 미지정 | 미지정 → 기본 **+2** | **v2가 AlleleSeq2 ASE 모드와 이미 같다.** §2.3 ② 정정 참조 |
| `twopassMode` | Basic (index 1개) | Basic (**index 2개에 독립 적용**) | §2.3 ③ 참조 |
| `sjdbOverhang` | 100 | 74 | 75bp read 기준 read_length−1. 타당 |
| `outSAMattributes` | All | `NH HI AS nM NM` | AS 포함 — `33`이 요구. 충족 |

### 2.3 실질적 문제 — 우선순위순

#### ① duplicate가 ASE 카운트에 그대로 들어간다 — 확실한 결함

`16_allelic_counts.py`의 `read_as()`는 `is_unmapped` / `is_secondary` / `is_supplementary` /
`mapping_quality`만 거르고 **`is_duplicate`를 확인하지 않는다.**
`31`이 MarkDuplicates를 `REMOVE_DUPLICATES=false`로 돌려 flag만 남기므로 DUP 리드가 독립 fragment로 카운트된다.

같은 파이프라인 안에서 `15_coverage_zarr.py`는 DUP를 제외하므로
**커버리지와 ASE 카운트의 duplicate 정책이 서로 다르다.** AlleleSeq2는 카운트 전에 duplicate를 제거한다
(`RM_DUPLICATE_READS=on` → `PIPELINE.mk`의 `rmdup.bam` 타깃).

영향: `totalCount` 부풀림 → 이항검정 과대 유의. 이는 docstring이 기록한 "read 단위 → fragment 단위"
수정(`is_ase` 0.45 문제)과 **같은 종류의 인플레이션**이며 아직 남아 있다.

```python
# 16_allelic_counts.py :: read_as()
if r.is_unmapped or r.is_secondary or r.is_supplementary:
    continue
# ← r.is_duplicate 체크 없음
```

#### ② `--scoreGenomicLengthLog2scale`을 0으로 두지 않았다

AlleleSeq2는 **ASE / ASB / ASCA / custom 전 모드에서 `--scoreGenomicLengthLog2scale 0.0`을 명시**한다
(STAR 기본값 −0.25). AS에 정렬 genomic length의 log2 비례 항이 섞이기 때문이다.

`31`은 이 옵션을 주지 않아 기본값 −0.25가 적용된다. HAP1/HAP2는 indel 때문에 같은 fragment의
genomic span이 다를 수 있고, 그 차이가 **염기 일치와 무관한 AS 차이**로 들어온다.
`33`은 `s1 == s2`를 tie로 버리고 1점 차이면 배정하므로 ±1 흔들림이 곧바로 배정/기각을 바꾼다.
AS 비교를 배정 규칙으로 쓰는 v2에서는 AlleleSeq2보다 오히려 더 critical하다.

> **정정 (소스 재확인).** `--scoreGenomicLengthLog2scale 0.0`이 4개 모드 전부에 있다는 것은 맞다.
> 그러나 초판이 함께 권고한 `--sjdbScore 0`은 **AlleleSeq2의 RNA 경로가 아니다.**
> `PIPELINE.mk`의 STAR 호출은 4개다:
>
> | 모드 | twopassMode | scoreGenomicLengthLog2scale | sjdbScore | scoreGap |
> |---|---|---|---|---|
> | `custom-params` | — | 0.0 | 기본 | — |
> | **`ASE-params`** ← v2 대응 | Basic | 0.0 | **기본(+2)** | — |
> | `ASB-params` | — | 0.0 | **0** | −100 |
> | `ASCA-params` | — | 0.0 | **0** | −100 |
>
> `sjdbScore 0`은 splicing이 없는 **DNA 모드(ASB/ASCA)** 전용이고 `--scoreGap -100`이 함께 붙는다.
> RNA-seq용 `ASE-params`는 기본값을 쓴다. 따라서 RNA 경로에 0을 넣는 것은 AlleleSeq2를 따르는 것이
> 아니라 그 DNA 모드를 RNA에 적용하는 것이다. **적용하지 않았다.**
> lift된 GTF junction의 비대칭 문제(아래 ③)는 실재하지만 그 해법은 sjdbScore가 아니라
> lift 정책 쪽이며, ⑨에서 다룬다.

#### ③ HAP1/HAP2에서 2-pass를 독립 실행한다

`--twopassMode Basic`이 두 index에서 따로 돌면서 pass-1에서 발견한 novel junction 집합이 **서로 달라진다.**
한쪽에만 삽입된 junction을 가로지르는 리드는 그쪽에서 AS가 크게 높아져 그 haplotype으로 배정된다 —
het와 무관한 신호다. AlleleSeq2는 index가 하나라 구조적으로 발생하지 않는다.

같은 이유로 `--sjdbScore`(기본 +2)도 문제다. lift된 GTF junction이 `lift_snapped`로 한쪽만 밀리면
+2가 비대칭으로 붙는다. AlleleSeq2가 ASB 모드에서 `--sjdbScore 0`을 주는 이유가 이것이다.
`13_lift_annotation.py`의 `gtf_lift_qc.json` `snap_fraction_of_output`을 HAP1/HAP2 간에 비교하면
이 비대칭의 크기를 실측할 수 있다.

#### ④ CNV 필터가 없다

AlleleSeq2는 WGS BAM에서 hetSNV ±1kb read depth를 계산해(`hetSNVs_rd.tab`)
**0.5 ≤ rd ≤ 1.5 밖이면 사이트를 버린다** (`filter_cnv_sites.py`).
v2에는 대응물이 전혀 없다. CNV 영역의 유전자는 구조적으로 allele-specific하게 보이므로
1000G 샘플에서는 확실한 위양성 소스다.

관련 AlleleSeq2 파일: `makePersonalGenome.mk`의 `hetSNVs_rd.tab` 타깃, `calculate_rd.py`,
`feed_sorted_check_names.py`, `filter_cnv_sites.py`.

#### ⑤ multimapping 보정이 없다

AlleleSeq2는 multimapper를 별도 BAM/mpileup으로 뽑아 `filter_sites_w_mmaps.py adjust`로
**약한 allele 쪽에 multimapper 카운트를 더해 불균형이 multimapping 탓인지 확인**한다.
v2는 MAPQ 255로 그냥 버린다. 버리는 것은 카운트에 보수적이지만 **편향을 막지는 못한다** —
강한 allele 쪽 자리만 uniquely mappable하면 불균형이 그대로 남는다.
`33`의 docstring이 주장하는 "재정렬로 개선되는 것은 배정량이 아니라 reference mapping bias"를
검증하려면 이 보정이 필요하다.

#### ⑥ beta-binomial(overdispersion 모델)이 없다

AlleleSeq2 / AlleleDB가 `alleledb_calcOverdispersion.R` + `alleledb_alleleseqBetabinomial.R`를 두는 이유는
**plain binomial이 overdispersion 때문에 ASE를 과대 호출하기 때문**이다.
`18_ase_qc.py`는 exact binomial + BH만 쓴다. split-half 상관으로 재현성을 경험적으로 보긴 하지만
통계 모델 차원의 보정은 아니다. `33`의 docstring이 기록한
"`is_ase` 0.45 (참조축 진리값 0.067, LCL 정상 범위 5–10%)"의 잔여분이 여기 있을 가능성이 높다.

#### ⑦ decoy / unplaced contig가 없다

`makePersonalGenome.mk`는 `ADDNL_SEQNS` (unplaced / unlocalized / EBV / decoy)를
**"siphon off reads" 목적으로** diploid index에 일부러 넣는다 (AS 계산에는 쓰지 않음).
v2는 `CHROMS=chr1-chr22`만 쓴다. chrX/Y/MT/unplaced 유래 리드가 갈 곳이 없어 autosome에 강제로 붙고,
그 오정렬은 두 hap에서 비대칭일 수 있다.

#### ⑧ 한쪽 hap에만 정렬된 fragment를 기본 포함한다

`only1` / `only2` (전체 배정의 10–19%)는 AS 비교가 아니라 정렬 성공/실패에서 나온 신호다.
`33`의 주석 자신이 "공통 mate가 없으면 정렬 성공/실패가 대립 증거로 둔갑한다"고 쓰면서
그 케이스(`n_nocommon`)는 버리지만, `only1`/`only2`는 기본 포함한다.
`.both_haps.tsv`와 `single_hap_sensitivity_delta`를 만들어 둔 것은 좋으나
**primary 산출물(`allelic_counts.tsv`)은 여전히 포함 버전**이다.

> **보충.** AlleleSeq2에는 이 구분 자체가 없다 — diploid 단일 index라 리드가 한 번만 정렬된다.
> 따라서 기본값을 뒤집을 근거를 AlleleSeq2에서 가져올 수 없다. 별도 판단으로 남긴다.

#### ⑨ lift 실패 정책이 정반대다 — 초판이 놓친 차이

AlleleSeq2는 hetSNV를 hap 좌표로 옮기지 못하면 `not_lifted?`로 표시하고
`filter_phase_warnings.py`가 **그 site를 버린다**.

v2의 `13_lift_annotation.py`는 정반대였다 — "가장 가까운 블록 경계로 스냅하고
`lift_snapped=1` 플래그를 남긴다 — 조용히 버리지 않는다"가 명시된 설계였다.

보존은 일반적으로 옳지만 ASE에서는 위험하다. **스냅은 HAP1/HAP2에서 비대칭으로 일어난다** —
두 haplotype의 결실 위치가 다르므로 같은 exon이 한쪽에서만 밀린다. 그러면 같은 fragment가
hap별로 다른 유전자에 귀속될 수 있고, 이는 het와 무관한 신호다. ③의 junction 비대칭과 같은
뿌리이며, 이쪽이 더 직접적이다.

#### ⑩ site 단위 품질 필터가 없다 — 초판이 놓친 차이

`filter_phase_warnings.py`는 lift 실패 외에 **개인 유전체가 그 자리에서 틀렸다는 신호**도 거른다.
`read_pileup.py`가 site마다 다음 warning을 붙이고, 하나라도 `.`이 아니면 site를 버린다:

- `N_other_alleles` — hap allele이 아닌 염기를 가진 read가 섞여 있다
- `hap_allele_not_largest_cnt` — hap allele이 최다 관측 염기가 아니다

v2는 fragment의 AS를 비교해 gene 단위로 세므로 **"이 자리의 개인 서열이 맞나"를 물어볼 지점이
구조적으로 없다.** 이식하려면 het SNV 위치를 hap 좌표로 옮기고 pileup을 뜨는 별도 경로가 필요하다.

### 2.4 문제없이 대응되는 항목

- **염색체 범위** — AlleleSeq2는 `filter_non-autosomal_chr.py`로 비-상염색체를 배제(`KEEP_CHR`로 X만 선택 허용).
  v2는 `CHROMS=chr1-chr22`로 동등.
- **indel 반영** — AlleleSeq2 vcf2diploid도 SNV+indel(+SV)을 개인 서열에 반영. v2의 `SEQUENCE_VARIANT_POLICY`
  (snv+indel)와 동등.
- **개인 좌표 주석 사용** — 양쪽 모두 참조 좌표 GTF를 그대로 쓰지 않고 chain으로 lift한다. `31`의 주석
  "참조 좌표 GTF를 그대로 주면 indel만큼 어긋난 junction이 삽입된다"는 AlleleSeq2와 같은 판단.
- **`sjdbOverhang`** — AlleleSeq2 100 vs v2 74. v2는 75bp read의 `read_length − 1`로 오히려 더 정확.
- **AS 태그 출력** — `--outSAMattributes`에 AS 포함. `33`의 요구 충족.

---

## 3. 조치 항목 — 2026-09-17 적용 결과

전부 이 디렉토리에서만 수정했다. 동결 원본 `mfm_data/process_code/v2/`는 건드리지 않았다.

### 적용함

| # | 조치 | 파일 | 결과 |
|---|---|---|---|
| 1 | duplicate 제외 (`--include-duplicates`로 복원 가능) | `16_allelic_counts.py` `read_as()` | 기본 제외. AlleleSeq2 `RM_DUPLICATE_READS=on`과 `15_coverage_zarr.py`에 정책을 맞췄다 |
| 2 | `--scoreGenomicLengthLog2scale 0` | `14_realign_rna_star.sh` | 적용. **`--sjdbScore 0`은 적용하지 않았다** — §2.3 ② 정정 |
| 5 | CNV 필터 | **`17_cnv_filter.py` 신규** + `18_ase_qc.py` | gene 단위 rd를 stage 22 개인 WGS 정렬에서 재고 `[0.5, 1.5]` 밖이면 `cnv_region`으로 제외 |
| 6 | beta-binomial | `18_ase_qc.py` | overdispersion rho를 AlleleSeq2와 같은 가중 SSE 히스토그램 적합으로 추정하고 `min(1, 2·pbetabinom)`로 양측 p를 만든다. **이쪽이 본선**이고 binomial은 비교용으로 함께 보고한다 |
| 7 | decoy / unplaced contig | `14_realign_rna_star.sh` | `ADDNL_SEQNS`(unplaced/unlocalized/EBV/decoy)를 **index에만** 넣는다. alt/HLA는 양쪽 모두 제외 |
| ⑤ | multimapping 보정 | `16_allelic_counts.py` | uniq(MAPQ≥255)와 mmap(0<MAPQ<255)을 서로소로 배정하고, 약한 쪽에 mmap을 몰아준 뒤에도 방향이 남는지 `.mmap_adjusted.tsv`에 적는다. 뒤집히면 `mmap_explains_imbalance` |
| ⑨ | lift 실패 시 drop | `13_lift_annotation.py` | `--on-unlifted drop`이 기본. `snap`으로 이전 동작 복원 가능 |

검증:
- beta-binomial pmf를 scipy `betabinom`과 대조 — 최대 오차 `7.6e-15`, rho→0에서 binomial과 일치.
- 합성 데이터로 rho 복원 — 참값 0.00 / 0.02 / 0.05 / 0.10 / 0.20을 격자 정밀도(0.01)로 모두 복원.
- CNV 필터를 합성 BAM으로 end-to-end — 정상 rd=1.0 통과, 3배 증폭 rd=3.0 제외, 결실 rd=0.2 제외.
  HAP1만 3배인 het 증폭도 제외(hap 하나라도 구간 밖이면 제외).
- 테스트 **13개 전부 통과** (`./tests/run_all.sh`).

### 적용하지 않음 — 근거

| # | 조치 | 이유 |
|---|---|---|
| 3 | HAP별 2-pass 분리 또는 SJ 집합 통합 | 설계 판단이 필요하다. AlleleSeq2는 단일 index라 대응물이 없어 "AlleleSeq2 기준"으로 정할 수 없다. 남은 가장 큰 미해결 항목이다 |
| 4 | `--require-both-haps`를 기본값으로 | **AlleleSeq2에는 only1/only2 구분 자체가 없다**(단일 index). 근거를 가져올 수 없어 기본값을 유지했다 |
| ⑩ | site 단위 품질 필터(`filter_phase_warnings.py`) | v2는 gene 단위 AS 비교라 site 판정 지점이 구조적으로 없다. het SNV 좌표 lift + pileup 경로를 새로 만들어야 한다 |
| — | 시뮬레이션 FDR(`FalsePos.py`, 500회) | BH를 유지했다. 난수 없이 재현되고 기존 컬럼과 이어진다. 과대 호출의 주원인인 overdispersion은 조치 6이 교정한다 |

### 스키마 변경 (기존 산출물과 구분)

| 대상 | 이전 | 현재 |
|---|---|---|
| `33` meta schema | `mfm_allelic_counts_meta_v2` | `mfm_allelic_counts_meta_v3` |
| `33` script_revision | `2026-09-16_fragment_split_half` | `2026-09-17_alleleseq2_dup_mmap` |
| `31` rna_realign schema | `mfm_v2_rna_realign_v2` | `mfm_v2_rna_realign_v3` |
| `30` lift QC schema | `mfm_gtf_lift_qc_v1` | `mfm_gtf_lift_qc_v2` |
| `35` QC schema | `mfm_ase_qc_v1` | `mfm_ase_qc_v2` |

`35`의 `provenance_ok`는 새 revision을 요구하므로 **기존 산출물은 `invalid_provenance`로 떨어진다.**
`31`의 skip 판정도 `v3`를 grep하므로 기존 BAM은 재사용되지 않고 재정렬된다. 의도된 동작이다.
`35`는 `include_duplicates=true`인 산출물도 `duplicates_included`로 막는다.

### 새 설정 (`config_v2.sh`)

```
RNA_INDEX_ADDNL_SEQNS=1        # decoy/unplaced/EBV를 STAR index에 포함
RNA_ADDNL_SEQNS_SOURCE=...     # 추출 원본 (Broad hs38DH)
RNA_ADDNL_SEQNS_FASTA=...      # 추출 결과 캐시 (최초 1회 생성)
RNA_ASE_INCLUDE_DUPLICATES=0   # 1이면 duplicate를 카운트에 포함
RNA_LIFT_ON_UNLIFTED=drop      # drop | snap
RNA_CNV_FILTER=1               # WGS 개인 정렬이 있을 때만 동작
RNA_CNV_MIN_RD=0.5
RNA_CNV_MAX_RD=1.5
```

⚠ `RNA_INDEX_ADDNL_SEQNS=1`은 STAR index 생성 비용을 늘린다(contig 약 2,550개 추가).
기존 index는 스키마 변경으로 어차피 재생성된다.

## 4. 참고

- AlleleSeq2 논문/도구: <http://alleleseq.gersteinlab.org/tools.html>
- 대조에 사용한 AlleleSeq2 파일:
  `PIPELINE.mk`, `makePersonalGenome.mk`, `pileup2counts.py`, `read_pileup.py`,
  `filter_cnv_sites.py`, `filter_phase_warnings.py`, `filter_sites_w_mmaps.py`,
  `filter_by_counts.py`, `filter_non-autosomal_chr.py`, `get_hetSNV_bed.py`,
  `FalsePos.py`, `alleledb_calcOverdispersion.R`, `alleledb_alleleseqBetabinomial.R`
- v2 내부 근거 문서 (`16_allelic_counts.py` docstring이 인용):
  `post-training/journal/2609/ase_truth_spread_260903.md` §9.4,
  `ase_validation_260820.md` §3, §4
