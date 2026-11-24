# integrated WGS/RNA 전처리 v2 인수인계

## 2026-09-16 현행 상태

이 절이 아래의 2026-08 개발 이력보다 우선한다. current cohort 및 외부 전달본은 `/NHNHOME/korea/mfm_data/process_code/v2/share_package`다.

| 구분 | 현재 값 |
|---|---:|
| 최종 고유 개인 | 3,212 |
| PacBio HiFi | 294 |
| PacBio CLR | 0 |
| Illumina 30x | 2,918 |
| strict GEUVADIS RNA BAM pair | 449 |
| WGS-only | 2,763 |

- PacBio 선택은 HPRC R2 232명, HGSVC3 60명, HGSVC2 1명, HPRC R1 1명이며 모두 HiFi다.
- HPRC R3는 release 변동성과 첫 production freeze의 재현성 때문에 제외했다.
- ONT 경로와 Clair3 self-call은 legacy/비교 실행으로 남아 있지만 current PacBio production cohort의 기본 자료원은 아니다.
- current PacBio 경로는 검증된 프로젝트별 public VCF를 우선하며 donor/reference/phasing 검증 실패 시 Clair3로 자동 fallback하지 않는다.
- manifest는 3,212명 모두의 WGS 형식·URL·bytes와 449명 RNA pair 정보를 기록한다. 실제 checksum 가용 여부는 source별 manifest 필드를 따른다.
- Stage 01은 aligned BAM/CRAM을 canonical BAM으로 표준화하고 raw HiFi uBAM/FASTQ는 minimap2 `map-hifi`로 정렬한다. 현재 production manifest에는 CLR/ONT가 없다.
- RNA branch는 STAR 재정렬, Picard mark-only, 전장/WGS-window coverage, validated-HAP ASE v5/v3, CNV, final validation과 cleanup까지 연결돼 있다.
- WGS-RNA pair의 중간 WGS BAM은 RNA/ASE와 post-training compact QC가 끝날 때까지 삭제하지 않는다.

## 위치

- v1 보존본: `/NHNHOME/korea/pre-training/jungjin/MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v1`
- v2 개발본: `/NHNHOME/korea/pre-training/jungjin/MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2`
- 실행 wrapper: `/NHNHOME/korea/pre-training/jungjin/MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2/scripts/run_integrated_preprocess.sh`
- 데이터 기준 root: `/NHNHOME/korea/mfm_data`
- 계획/설계 기록: `/NHNHOME/korea/pre-training/jungjin/codex_archive/plan/260825`, `/NHNHOME/korea/pre-training/jungjin/codex_archive/plan/260826`
- 코드 백업: `/NHNHOME/korea/pre-training/jungjin/codex_archive/code_backup/260825`, `/NHNHOME/korea/pre-training/jungjin/codex_archive/code_backup/260826`

## v2 핵심 정책

1. long-read production 기본값은 manifest의 검증된 donor-matched public phased VCF다. Clair3 self-call은 VCF가 없고 명시적으로 승인한 비교 실행에서만 사용하며 검증 실패의 자동 fallback이 아니다.
2. 공개 VCF 분기는 VCF 생성에 사용된 정확한 reference를 CRAM/BAM 검증과 consensus에 함께 사용한다. contig·M5·REF allele 불일치는 자동 보정하지 않고 실패 처리한다.
3. short-read Illumina는 1000G Illumina VCF 기반 FASTA 생성과 WhatsHap haplotagging을 허용한다.
4. long-read haplotype 정보는 Clair3 VCF 또는 검증된 public phased VCF를 WhatsHap `haplotag` 입력으로 사용해 BAM에 `HP`/`PS`를 부여한다. 직접 BAM-only clustering은 fallback/비교 실험용이다.
5. long-read에서 1000G Illumina panel VCF로 FASTA를 만들거나 haplotagging하는 경로는 제거했다.
6. final WGS compact는 HAP별 personal FASTA 좌표와 personal-realigned readset을 기준으로 만든다.
7. final compact에는 per-base depth를 `uint32`로 함께 저장한다.
8. 기본 `HAP_READ_MODE=all_to_both`는 read를 양 HAP에 재정렬해 coverage를 보존한다. `hp_routed_with_unknown`과 `hap_specific`은 목적에 따라 명시적으로 선택하며 mode는 layout/provenance에 기록한다.
9. GEUVADIS RNA는 SNV+indel personal FASTA에 STAR로 재정렬하고, duplicate mark·전장 및 WGS-window coverage·validated-HAP ASE·CNV·metadata·최종 검증을 수행한다.

## 단계 구조

### 01 cohort common VCF

원본 1000G chromosome panel VCF에서 cohort sample만 subset한다.

```bash
bash MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2/scripts/run_integrated_preprocess.sh \
  extract-cohort-common-vcf \
  --cohort <cohort_id> \
  --sample-list <samples.txt> \
  --output-dir <work/common_vcf> \
  --chroms chr1-chr22 \
  --panel-base <1000G_panel_dir_or_url> \
  --panel-file-template '1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz' \
  --filter-status PASS,. \
  --min-ac 1 \
  --jobs 6 \
  --threads-per-chrom 2
```

### 02 per-sample VCF

sample별 VCF를 만든다.

- `<sample>.seq.vcf.gz`: personal FASTA consensus용. v2 기본 `snv_indel`.
- `<sample>.tag.vcf.gz`: short-read WhatsHap haplotagging용. long-read production에서는 사용 금지.

```bash
bash MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2/scripts/run_integrated_preprocess.sh \
  extract-per-sample-vcfs \
  --cohort <cohort_id> \
  --common-vcf-dir <work/common_vcf> \
  --output-dir <work/persample_vcf> \
  --sample-list <samples.txt> \
  --chroms chr1-chr22 \
  --sequence-variant-policy snv_indel \
  --jobs 8
```

### 03 personal FASTA

HAP1/HAP2 FASTA와 chain을 생성한다.

```bash
bash MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2/scripts/run_integrated_preprocess.sh \
  create-personal-fastas \
  --cohort <cohort_id> \
  --seq-vcf-dir <work/persample_vcf> \
  --output-root <final/pre_training/{read_type}/{project_id}> \
  --reference-fasta <ucsc_hg38_autosomes.fa> \
  --sample-list <samples.txt> \
  --chroms chr1-chr22 \
  --sequence-variant-policy snv_indel \
  --jobs 8
```

출력:

```text
<final-root>/<sample>/individual_fasta/HAP1/chr*.fa, chr*.fa.fai, chr*.chain
<final-root>/<sample>/individual_fasta/HAP2/chr*.fa, chr*.fa.fai, chr*.chain
```

### 04 short-read haplotagging

short-read Illumina만 기본 허용한다.

```bash
bash MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2/scripts/run_integrated_preprocess.sh \
  create-haplotagged-alignments \
  --cohort <cohort_id> \
  --raw-root <raw_root> \
  --tag-vcf-dir <work/persample_vcf> \
  --output-root <final-root> \
  --reference-fasta <vienna_or_decode_ref.fa> \
  --platform Illumina \
  --sample-list <samples.txt>
```

long-read에서 이 명령을 직접 사용하면 기본적으로 실패한다. stage 10에서 검증된 public VCF에 한해 wrapper가 `--allow-longread-vcf-haplotag`을 명시해 호출한다. 기본 Clair3 source는 `04_call_longread_clair3.sh`에서 처리한다.

### 04L long-read Clair3 self-calling + WhatsHap haplotagging

ONT/HiFi long-read에서는 Illumina/statistical VCF로 직접 haplotagging하지 않는다. 대신 같은 long-read BAM/CRAM에서 Clair3로 variant를 call하고, Clair3가 만든 phased VCF를 WhatsHap `haplotag` 입력으로 사용해 원본 BAM/CRAM에 `HP`/`PS` tag를 부여한다.

이 단계에서 WhatsHap을 사용하는 이유는 명확하다. WhatsHap `haplotag`는 이미 phased VCF가 있을 때 각 read가 어느 haplotype에 더 잘 맞는지 판정하고 BAM에 `HP:i:1/2`와 `PS:i:<phase_set>` tag를 붙이는 표준적인 역할을 한다. 단, 이 VCF는 Illumina VCF가 아니라 Clair3가 long-read에서 새로 call/phasing한 VCF여야 한다.

```bash
bash MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2/scripts/run_integrated_preprocess.sh \
  call-longread-clair3 \
  --cohort <cohort_id> \
  --output-root <final-root> \
  --reference-fasta <cram_decode_reference.fa> \
  --model-path <clair3_model_dir> \
  --platform ont \
  --raw-root <raw_root> \
  --sample-list <samples.txt> \
  --alignment-pattern '{sample}/wgs/{sample}.hg38.cram' \
  --contigs chr1-chr22 \
  --threads 32 \
  --output-threads 8 \
  --jobs 1
```

단일 sample/chr22 smoke 예시:

```bash
bash MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2/scripts/run_integrated_preprocess.sh \
  call-longread-clair3 \
  --cohort smoke \
  --output-root /tmp/clair3_v2_smoke \
  --reference-fasta /path/to/vienna_or_decode_ref.fa \
  --model-path /path/to/clair3_model \
  --platform ont \
  --alignment /path/to/HG00100.hg38.cram \
  --sample HG00100 \
  --contigs chr22 \
  --threads 32 \
  --output-threads 8
```

출력:

```text
<final-root>/<sample>/longread_calling/clair3/
  merge_output.vcf.gz
  ... Clair3 산출물 ...

<final-root>/<sample>/haplotagged/
  <sample>.clair3_longread.haplotagged.bam
  <sample>.clair3_longread.haplotagged.bam.bai
  <sample>.clair3_longread.taglist.tsv.gz

<final-root>/<sample>/qc/
  clair3_longread_haplotag.json
```

검증 기준:

- `merge_output.vcf.gz` 존재 및 index 생성 가능 여부
- haplotagged BAM 존재 및 index 존재 여부
- inspected read 중 `HP`/`PS` tag 비율
- 기존 직접 BAM-only 결과 대비 unknown/unphased 감소 여부
- Illumina VCF와 phase concordance QC. Illumina VCF는 입력이 아니라 비교 기준이다.

### 04X long-read BAM-only phasing fallback / read clustering

이 단계는 production 기본 경로가 아니다. Clair3/WhatsHap이 설치 문제, 모델 문제, 특정 실험 조건 문제로 막힐 때 비교용 또는 fallback으로만 사용한다.

기존 직접 구현된 방식:

현재 구현된 상태:

- BAM/CRAM에서 C-backed `pysam.count_coverage` chunk 방식으로 biallelic SNV-like marker를 찾는다.
- read별 marker allele observation을 만든다. 기본 baseQ filter는 0이며, baseQ가 유효한 BAM에서만 높인다.
- marker gap/marker 수 기준으로 phase block 후보를 나눈다.
- allele pattern으로 read를 `cluster=1`, `cluster=2`, `cluster=0(unknown)`으로 분류한다.
- `cluster=1/2`는 아직 HAP1/HAP2와 대응되지 않는다.
- HP/PS-tagged BAM/CRAM writer는 아직 없다.

```bash
bash MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2/scripts/run_integrated_preprocess.sh \
  phase-longread-bamonly \
  --sample-id <sample> \
  --alignment <raw_longread.hg38.cram_or_bam> \
  --output-dir <final-root>/<sample>/bamonly_phase \
  --reference-fasta <cram_decode_reference.fa> \
  --platform ont \
  --regions chr22:20000000-21000000
```

산출물:

```text
<sample>.bamonly_phase.markers.tsv.gz
<sample>.bamonly_phase.read_clusters.tsv.gz
<sample>.bamonly_phase.qc.json
```

주의:

- 이 결과를 바로 `HP=1/2`로 간주하면 안 된다.
- 다음 구현에서 cluster-to-HAP orientation이 필요하다.
- orientation 전 long-read personal realignment 실험은 `--read-routing-policy all_to_both`로 수행할 수 있다.

### 05 window index

HAP1/HAP2 FASTA에서 1M window를 만든다. 시작점은 seed/sample/haplotype/contig 기반 deterministic random offset이며 기본 offset 범위는 `[0, 1M)`이다.

```bash
bash MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2/scripts/run_integrated_preprocess.sh \
  create-window-index \
  --cohort <cohort_id> \
  --output-root <final-root> \
  --sample-list <samples.txt> \
  --chroms chr1-chr22 \
  --sequence-variant-policy snv_indel \
  --window-size 1048576 \
  --seed 42 \
  --jobs 8
```

### 06 personal realignment

WGS read를 HAP1/HAP2 FASTA에 다시 정렬한다. ONT는 `map-ont`, HiFi/CCS는 `map-hifi`, Illumina는 `sr` preset을 자동 선택할 수 있다.

Clair3/WhatsHap HP/PS BAM을 아직 만들지 않은 상태에서의 long-read fallback:

```bash
bash MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2/scripts/run_integrated_preprocess.sh \
  create-personal-realigned-alignments \
  --cohort <cohort_id> \
  --output-root <final-root> \
  --sample <sample> \
  --chroms chr1-chr22 \
  --sequence-variant-policy snv_indel \
  --dataset-role pre_training \
  --alignment <raw_longread.hg38.cram_or_bam> \
  --reference-fasta <cram_decode_reference.fa> \
  --platform ont \
  --preset auto \
  --read-routing-policy all_to_both \
  --unknown-read-policy include_both \
  --format cram \
  --threads 24 \
  --sort-threads 8
```

Clair3/WhatsHap으로 HP/PS BAM이 생기면 `--read-routing-policy hp_routed`로 전환한다.

### 07 final compact

모델 입력용 zarr를 만든다. v2 기본 layout은 `full_shared_hap_4bit_acgtnp_depth`다.

주요 저장 데이터:

- `hap_sequence.zarr`: HAP sequence compact 4-bit A/C/G/T/N/P
- HAP별 `readset.zarr`: read sequence, CIGAR/mapping/clipping/strand/HP/PS 등 read metadata
- `depth.zarr`: window별 per-base depth, `uint32`
- window metadata: HAP coordinate, source reference info, offset 정보

검증:

```bash
bash MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2/scripts/run_integrated_preprocess.sh \
  validate-compact-v2 \
  --compact-contig-dir <final-root>/<sample>/model_input/full_shared_hap_4bit_acgtnp_depth/chr22 \
  --read-coordinate-mode personal_realigned \
  --validation-level standard \
  --output-json <qc.json>
```

## 현재 남은 구현 과제

1. production 규모 long-read Clair3 실행 시간/HP·PS tag 비율 측정
   - smoke에서는 작은 slice와 fallback을 사용했으므로, 실제 sample chr22 전체 또는 autosome 전체에서 phased VCF 생성 여부를 확인해야 한다.
2. Clair3-derived phased VCF vs Illumina VCF QC
   - Illumina VCF는 long-read tagging 입력이 아니라 phase concordance QC 기준으로만 사용한다.
3. BAM-only fallback 비교 QC
   - 기존 직접 구현 방식의 unknown 비율과 Clair3/WhatsHap 결과의 untagged 비율을 비교한다.
4. GEUVADIS RNA branch production smoke
   - 기존 process_code RNA branch는 보존했다. WGS compact 내부 RNA coverage 통합은 후속 구현이다.
5. 대규모 production smoke
   - 샘플 1개 chr22 전체
   - 샘플 1개 autosome 전체
   - cohort sample list 병렬 실행

## 검증 완료 항목

- v2 Python module `py_compile`
- `phase-longread-bamonly --help`
- synthetic read allele pattern clustering smoke
- dry-run contract smoke

## process_code형 root wrapper 구조

`/NHNHOME/korea/mfm_data/process_code/v2`와 병합하기 쉽도록 root-level 단계별 wrapper를 추가했다. 기존 Python backend는 `src/integrate_wgs_rna`에 유지하고, root `.sh` 파일은 해당 CLI를 호출한다. shell은 단계 orchestration, Python package는 실제 구현으로 역할을 분리한다.

```text
00_environment.sh          # conda-only 환경 설치/검증
config_v2.sh               # 공통 설정/manifest parser
03_prepare_vcf.sh          # common/per-sample VCF 준비
04_call_longread_clair3.sh # long-read Clair3 self-calling + WhatsHap haplotag
06_haplotag.sh             # long-read 산출물 확인 / short-read WhatsHap haplotag
09_personal_fasta.sh       # HAP1/HAP2 personal FASTA + chain
10_realign_wgs.sh          # WGS read personal realignment
11_windows.sh              # HAP별 1M window index
12_compact_zarr.sh         # final compact model input zarr
run_wgs_sample.sh          # sample 1개 WGS 전체 실행
runner_local.sh            # local multi-sample runner
slurm_wgs.sbatch           # SLURM WGS array template
GUIDE.md                   # 실행 가이드
examples/sample_manifest.tsv
```

process_code/v2 병합본에서는 기존 STAR personal realignment RNA branch를 보존한다. 자세한 실행법은 `GUIDE.md`를 기준으로 한다.



### 00_environment.sh 주요 모드

```bash
./00_environment.sh verify         # core + Clair3 + reference 검증
./00_environment.sh verify-core    # jj_mfm/core tool 검증
./00_environment.sh verify-clair3  # Clair3 전용 env 검증
./00_environment.sh verify-rna     # RNA tool 검증
./00_environment.sh install-core   # jj_mfm core package 설치/보강
./00_environment.sh install-clair3 # Clair3 전용 env 설치/보강
./00_environment.sh install-rna    # STAR/bedtools/pyBigWig 등 RNA tool 설치/보강
./00_environment.sh install-all    # 전체 설치/보강
```

현재 WGS core/Clair3 검증과 WGS smoke-wgs는 통과한다. RNA branch는 기존 process_code 구현을 보존했으며, `jj_mfm`의 STAR/Picard/OpenJDK/bedtools/samtools 구성으로 `verify-rna`를 통과한다. RNA production sample별 personal alignment·ASE QC는 별도 수행이 필요하다.


## 대표 bash `run_pipeline.sh`

단계별 스크립트를 직접 호출하지 않고 다음 entrypoint 하나로 실행할 수 있다.

```bash
./run_pipeline.sh env verify
./run_pipeline.sh wgs <sample>
./run_pipeline.sh wgs-list <sample_list>
./run_pipeline.sh slurm-wgs <sample_list> '1-20%4'
./run_pipeline.sh smoke
./run_pipeline.sh smoke-wgs <sample>
```

`smoke`는 빠른 계약 검증이고, 실제 작은 전처리는 `smoke-wgs`에서 수행한다.

## 2026-08-26 process_code/v2 병합 결과

병합 위치는 `/NHNHOME/korea/mfm_data/process_code/v2`다. integrated v2 WGS/core backend를 병합했고, 기존 process_code RNA branch는 삭제하지 않고 유지했다. 대표 실행은 `run_pipeline.sh` 하나로 통합했다.

검증 결과:

```text
./run_pipeline.sh smoke              PASS
./run_pipeline.sh smoke-wgs HG00096  PASS
```

`smoke-wgs HG00096`는 `chr22:20,000,000-21,000,000` slice를 사용해 다음 전체 경로를 통과했다.

```text
10 per-sample VCF
20 Clair3 long-read call + WhatsHap haplotag
21 haplotagged BAM 확인
40 HAP1/HAP2 personal FASTA
41 HAP별 1M window index
22 HAP personal realignment
70 compact zarr
```

smoke output 핵심값:

```text
compact layout: full_shared_hap_4bit_acgtnp_depth
readset layout: hap_split
read count: 3,536
HAP1 windows/segments: 28 / 2,297
HAP2 windows/segments: 27 / 1,239
depth dtype: uint32
smoke tmp size: 약 350 MB
```

주의: smoke에서는 작은 BED/sliced BAM 조건 때문에 Clair3 `pileup.vcf` fallback과 unphased compact 허용 옵션을 켰다. production 기본값에서는 이 예외가 꺼져 있다.

## RNA branch 현재 실행 조건

WGS 경로는 `smoke`와 `smoke-wgs`를 통과했다. RNA branch는 기존 process_code STAR/coverage 구현을 보존했고, 현재 `jj_mfm` 환경에는 STAR/Picard/OpenJDK/bedtools/samtools가 설치되어 `./run_pipeline.sh env verify-rna`를 통과한다. RNA production sample별 personal alignment·ASE QC는 별도 검증한다.

```text
missing RNA executable: /home/aibio_1/miniconda3/envs/jj_mfm/bin/bedtools
```

RNA 전처리까지 한 번에 실행하려면 먼저 다음 중 하나가 필요하다.

```bash
cd /NHNHOME/korea/mfm_data/process_code/v2
./run_pipeline.sh env install-rna
# 또는 별도 RNA env를 쓰는 경우 STAR_BIN, BEDTOOLS_BIN, SAMTOOLS_BIN, PYTHON_BIN 경로를 명시
```

이 제한은 WGS 전처리 성공 여부와 별개다.



## 2026-08-26 HG00096 전체 autosome WGS 테스트

`/NHNHOME/korea/mfm_data/process_code/v2` 병합본으로 HG00096 ONT/GEUVADIS WGS CRAM을 사용해 `chr1-chr22` 전체 autosome 전처리를 실행했다.

### 테스트 입력

```text
sample: HG00096
readset: /NHNHOME/korea/mfm_data/igsr_ont_geuvadis/HG00096/wgs/HG00096.hg38.cram
CRAM reference: /NHNHOME/korea/pre-training/data/reference/1KG_ONT_VIENNA/1KG_ONT_VIENNA_hg38.fa
consensus reference: /NHNHOME/korea/mfm_data/reference/hg38_autosomes.fa
common VCF: /NHNHOME/korea/mfm_data/preprocessed_data/work/260813/igsr_ont_geuvadis/common_vcf
output root: /NHNHOME/korea/mfm_data/preprocessed_data/test_runs/260826_full_HG00096
chromosomes: chr1-chr22
```

### 실행 결과

전체 WGS 경로 `10 → 20 → 21 → 40 → 41 → 22 → 70`이 통과했다. 최종 compact validation은 `PASS 22/22`다.

```text
sample final root: /NHNHOME/korea/mfm_data/preprocessed_data/test_runs/260826_full_HG00096/final/pre_training/long_read/full_test_ont_geuvadis/HG00096
compact root: /NHNHOME/korea/mfm_data/preprocessed_data/test_runs/260826_full_HG00096/final/pre_training/long_read/full_test_ont_geuvadis/HG00096/compact/full_shared_hap_4bit_acgtnp_depth
sample final size: 약 81G
compact size: 약 20G
```

주요 compact 통계:

```text
contigs: 22
shared_reads: 10,495,848
HAP1_windows: 2,311
HAP2_windows: 2,309
HAP1_segments: 5,283,825
HAP2_segments: 5,272,748
HAP1_mean_depth: 13.6878
HAP2_mean_depth: 13.6540
HAP1 phased reads: 3,204,539
HAP2 phased reads: 3,194,240
```

### 테스트 중 수정된 production 코드

1. Clair3 subprocess PATH 보강
   - 파일: `src/integrate_wgs_rna/longread_clair3.py`
   - 이유: Clair3 내부 sort/compress 단계에서 `bgzip`, `tabix`를 찾지 못하면 `merge_output.vcf.gz`가 생성되지 않을 수 있다.
   - 처리: absolute executable parent directory를 subprocess `PATH`에 추가한다.

2. Clair3 overwrite 옵션 노출
   - 파일: `04_call_longread_clair3.sh`
   - 이유: 실패한 Clair3 작업을 같은 output directory에서 안전하게 재실행할 수 있어야 한다.
   - 처리: `CLAIR3_OVERWRITE=1`이면 `--overwrite`를 전달한다.

3. depth padding 정합성 수정
   - 파일: `src/integrate_wgs_rna/direct_compact_writer.py`
   - 이유: 마지막 partial window는 sequence가 `P`로 padding되어 저장되므로 depth 배열도 저장 sequence 길이와 같아야 한다.
   - 처리: depth storage length는 `len(window.sequence)` 기준으로 만들고, 실제 coverage는 `actual_length` 영역에만 누적한다. padding 영역 depth는 0이다.

동일 수정은 integrated 개발본 `MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA_v2`에도 반영되어 있다.

### 아직 삭제하지 않은 테스트 산출물

- 전체 샘플 테스트 output: `/NHNHOME/korea/mfm_data/preprocessed_data/test_runs/260826_full_HG00096`
- 이전에 중단한 chr22 smoke output: `/NHNHOME/korea/mfm_data/preprocessed_data/test_runs/260826_chr22_HG00096`

두 경로 모두 production 경로가 아니며, 삭제는 별도 확인 후 수행한다.


## 2026-08-26 long-read de novo assembly 개발 경로 추가

정책 변경 방향:

- short-read Illumina 샘플은 기존처럼 Illumina phased VCF로 HAP1/HAP2 personal FASTA를 만들고, short-read BAM/CRAM은 Illumina VCF로 WhatsHap haplotagging한다.
- long-read ONT/HiFi 샘플은 read haplotagging에는 Clair3 self-called VCF를 사용한다. 이 부분은 v2에 이미 반영되어 있다.
- 새로 바뀌는 핵심은 long-read sequence FASTA다. 목표는 short-read VCF consensus가 아니라 long-read BAM/CRAM 기반 de novo assembly FASTA를 사용하는 것이다.

추가된 코드:

- `07_denovo_assembly.sh`
- `08_anchor_assembly.sh`
- `src/integrate_wgs_rna/denovo_assembly.py`
- `src/integrate_wgs_rna/assembly_anchor.py`
- CLI command: `run-denovo-assembly`, `anchor-denovo-assembly`
- env mode: `./run_pipeline.sh env install-denovo`, `./run_pipeline.sh env verify-denovo`

hifiasm 설치:

```text
path: /home/aibio_1/miniconda3/envs/jj_mfm/bin/hifiasm
version: 0.25.0-r726
```

HG00096 chr22 benchmark:

```text
25_denovo_assembly elapsed: 82.6 sec
hifiasm real time: 78.3 sec
hifiasm peak RSS: 17.1 GB
FASTQ reads: 132,815
cleanup 전 FASTQ: 334.8 MB
cleanup 전 hifiasm .bin: 약 915 MB
cleanup 후 output root: 209 MB
26_anchor_assembly elapsed: 175.2 sec
```

HAP assembly 결과:

```text
hap1 p_ctg: 70 contigs, 14.19 Mb, anchored fraction 0.9829
hap2 p_ctg: 26 contigs, 6.22 Mb, anchored fraction 0.9794
primary p_ctg: 85 contigs, 16.71 Mb, anchored fraction 0.9767
```

해석:

- chr22 subset에서도 hifiasm은 hap1/hap2 p_ctg GFA를 생성했다.
- 하지만 길이가 chr22 전체를 대표하기에는 부족하므로, 이 subset 결과를 바로 production HAP FASTA로 쓰면 안 된다.
- 현재 구현은 assembly benchmark와 anchoring QC까지다. chromosome-scale scaffolding 또는 contig-coordinate model input 정책은 아직 결정해야 한다.

운영 기본값:

```bash
SEQUENCE_SOURCE=vcf_consensus
```

기존 one-click WGS pipeline을 깨지 않기 위해 기본값은 기존 경로로 유지했다. de novo 경로는 `SEQUENCE_SOURCE=denovo` 또는 `SEQUENCE_SOURCE=auto`를 명시할 때만 실행된다.

## long-read VCF source 정책

기본 WGS FASTA 생성 방식은 VCF-consensus다. long-read 기본값은 같은 BAM/CRAM의 Clair3 self-call이며, 명시적 manifest에서는 검증된 donor-matched 공개 phased VCF도 사용할 수 있다. 검증 실패 시 Clair3로 자동 fallback하지 않는다.

플랫폼별 VCF source는 다음과 같다.

| 플랫폼 | FASTA 생성 VCF | haplotagging VCF/BAM | 비고 |
|---|---|---|---|
| Illumina short-read | 1000G Illumina phased panel에서 추출한 `<sample>.seq.vcf.gz` | 1000G panel 기반 `<sample>.tag.vcf.gz` + WhatsHap | 3,202명 panel 일관성 유지 |
| ONT/HiFi long-read 기본 | 같은 long-read BAM/CRAM의 Clair3 VCF | Clair3 VCF + WhatsHap | `platform_matched` 기본값 |
| ONT/HiFi long-read 공개 VCF | 검증·정규화된 donor-matched phased SNP/indel VCF | 같은 공개 VCF + WhatsHap | `manifest` 명시, 정확한 VCF reference와 동일 reference alignment 필수 |
| long-read de novo | 사용자가 `SEQUENCE_SOURCE=denovo`를 명시한 경우만 stage 25/26 실행 | 별도 | 현재 production FASTA로 자동 승격하지 않음 |

long-read 기본 실행 흐름은 다음과 같다.

1. `03_prepare_vcf.sh`: Clair3 source는 skip하고, public source는 donor·contig·reference·phasing을 검증해 seq/tag VCF를 만든다.
2. `04_call_longread_clair3.sh`, `05_clair3_qc.sh`: Clair3 source에서만 실행한다.
3. `06_haplotag.sh`: 선택된 source의 phased tag VCF로 WhatsHap HP/PS tag를 만든다.
4. `09_personal_fasta.sh`: 선택된 seq VCF와 정확한 reference로 HAP1/HAP2 FASTA를 만든다.

주요 옵션:

- `VCF_SOURCE_POLICY=platform_matched|manifest`: 기본은 `platform_matched`, public 사용 시 `manifest`다.
- `LONGREAD_SEQUENCE_VCF_SOURCE=auto|clair3|public`: manifest 없는 단일-sample 실행용 override다.
- `LONGREAD_HAPLOTAG_VCF_SOURCE=auto|clair3|public`: sequence source와 동일해야 하며 상충하면 실패한다.
- `LONGREAD_SEQ_VCF_EXPORT_MODE=copy|symlink`: Clair3 VCF를 `<sample>.seq.vcf.gz` 위치로 내보내는 방식. 기본은 `copy`.
- `LONGREAD_ALLOW_UNPHASED_SEQUENCE_VCF=1`: 진단용. phased heterozygous GT가 없는 VCF로도 FASTA 생성을 허용한다. 기본은 비허용.

주의사항:

- Clair3 calling에는 `REF_DECODE`를 사용하고, personal FASTA consensus에는 `REF_SEQUENCE`를 사용한다. 두 reference가 다를 수 있으므로 stage 40의 `bcftools norm --check-ref x`가 mismatch variant를 제거하며, `dropped_ref_mismatch_count`를 QC로 확인해야 한다.
- Clair3 VCF에 phased heterozygous GT가 전혀 없으면 HAP1/HAP2 FASTA consensus가 안전하지 않으므로 기본적으로 실패 처리한다.

## 260827 정책 변경: read routing 기본값 all_to_both

WGS read를 HAP1/HAP2 personal FASTA에 재정렬할 때 기본값을 다음으로 변경했다.

```bash
READ_ROUTING_POLICY=all_to_both
UNKNOWN_READ_POLICY=include_both
```

의미:

- 모든 mapped primary read를 HAP1 FASTA와 HAP2 FASTA 양쪽에 각각 minimap2로 재정렬한다.
- untagged/HP0 read도 양쪽에 포함된다.
- unmapped read도 기존처럼 양쪽에 포함되며 `XU:i:1`로 표시된다.
- HP/PS tag는 read metadata로 유지되므로, 모델 또는 QC에서 원래 haplotype assignment를 확인할 수 있다.

주의:

- HAP1/HAP2 readset에 같은 read가 중복으로 존재할 수 있다.
- HAP별 depth는 biological haplotype-specific depth가 아니라 `per-haplotype realignment evidence depth`로 해석해야 한다.
- read 누락을 줄이고 모델링 유연성을 확보하기 위한 기본값이다.
- 생물학적으로 엄격한 HAP별 read 분리가 필요한 ablation/QC에서는 `HAP_READ_MODE=hap_specific`을 명시한다.

## 260828 고병렬 운영 정책 추가

대량 전처리는 `HIGH_PARALLEL_POLICY.md`를 따른다. 핵심은 다음과 같다.

- GPU 기본 미사용
- `LOCAL_JOBS=all|max` 지원
- sample-level 병렬 우선
- sample당 `CLAIR3_THREADS=12`, `WGS_REALIGN_THREADS=8`, `WGS_REALIGN_JOBS=2`, `COMPACT_JOBS=16` 권장
- 최종 compact 생성 후 중간 산출물 삭제 유지

현재 random20 실행 스크립트(`/NHNHOME/korea/pre-training/jungjin/codex_archive/plan/260827/run_random20_process_v2.sh`)도 다음 실행부터 `LOCAL_JOBS_LONG=32`를 기본값으로 사용한다. 이미 실행 중인 프로세스에는 소급 적용되지 않으므로, 현재 실행 배치는 기존 설정으로 계속 진행한다.
## 32개 배치 운영 및 sample_count 옵션

실제 운영 기본 sample-level 병렬도는 32개이다. `runner_local.sh`의 기본 `LOCAL_JOBS`가 32로 설정되어 있으며, 메모리보다 CPU와 Lustre I/O가 먼저 병목이 되므로 무제한 동시 실행은 실험용으로만 사용한다.

```bash
# 목록 전체를 최대 32개씩 동시에 실행
./run_pipeline.sh wgs-list /path/to/samples.txt

# 목록의 유효한 샘플 중 앞에서부터 10개만 실행
./run_pipeline.sh wgs-list /path/to/samples.txt 10

# RNA도 동일한 방식
./run_pipeline.sh rna-list /path/to/rna_samples.txt 10

# 동시 실행 수를 명시적으로 변경
LOCAL_JOBS=16 ./run_pipeline.sh wgs-list /path/to/samples.txt 10

# 실험용: 목록 전체를 동시에 시작
LOCAL_JOBS=all ./run_pipeline.sh wgs-list /path/to/samples.txt
```

`sample_count`는 주석(`#`)과 빈 줄을 제외한 목록 순서 기준 앞 N개를 의미한다. 무작위 표본이 필요하면 실행 전에 별도 목록 파일을 만들어 사용한다. 현재 실행 중인 작업에는 소급 적용되지 않고, 다음 실행부터 적용된다.



## post-training 처리 정책(최종)

post-training GEUVADIS 샘플도 pre-training과 동일하게 WGS 전처리(Clair3/WhatsHap, 개인 FASTA, HAP1/HAP2 realignment, window, compact zarr)를 수행한다. 이후 RNA만 추가로 STAR realignment와 1bp coverage/allelic count를 생성한다. WGS cleanup은 RNA 산출물이 완성될 때까지 자동 유예한다.


<!-- MFM_HIFI_QC_UPDATE_260828 -->
## HiFi/ONT Clair3 모델 선택 및 QC 정책

- long-read ONT/HiFi 기본값은 `VCF_SOURCE_POLICY=platform_matched`와 Clair3 self-call이다. 공개 phased VCF는 `VCF_SOURCE_POLICY=manifest`를 명시하고 검증을 통과한 sample에서만 사용한다.
- Clair3 모델은 샘플 플랫폼으로 자동 선택함. ONT는 `CLAIR3_MODEL_PATH_ONT` 기본값인 `$CLAIR3_MODELS_DIR/ont`, HiFi/CCS는 `CLAIR3_MODEL_PATH_HIFI` 기본값인 `$CLAIR3_MODELS_DIR/hifi`를 사용함.
- `CLAIR3_MODEL_PATH`를 직접 지정하면 기존 호환성을 위해 모든 플랫폼 공통 override로 처리함. 플랫폼별로 바꾸려면 `CLAIR3_MODEL_PATH_ONT`, `CLAIR3_MODEL_PATH_HIFI`를 지정함.
- WGS long-read pipeline은 `04_call_longread_clair3.sh` 직후 `05_clair3_qc.sh`를 실행함. QC는 `<final_root>/<sample>/qc/clair3_output_qc.json` 및 `.tsv`에 저장함.
- QC 항목은 모델-플랫폼 매칭, Clair3 VCF total/PASS/SNV/indel/other/multiallelic, GT/phased het fraction, haplotagged BAM의 HP/PS tag fraction, `min_mq`, `min_coverage`, `qual` 설정값임.
- QC만 재실행하려면 `./run_pipeline.sh clair3-qc <sample>` 또는 `LOCAL_JOBS=16 ./run_pipeline.sh clair3-qc-list <sample_list> [count]`를 사용함.

## 260917 HAP read mode

Stage 22는 HAP_READ_MODE로 read routing의 목적을 명시한다. 이 변수는 기존 READ_ROUTING_POLICY와 UNKNOWN_READ_POLICY보다 우선한다.

~~~bash
# 기본: coverage/evidence 보존
HAP_READ_MODE=all_to_both ./10_realign_wgs.sh <sample>

# HP0/missing을 양쪽에 보존하되 HP-tagged read는 해당 HAP에 보존
HAP_READ_MODE=hp_routed_with_unknown ./10_realign_wgs.sh <sample>

# biological HAP-specific 모델 입력용 별도 재실행
HAP_READ_MODE=hap_specific ./10_realign_wgs.sh <sample>
~~~

- all_to_both: 모든 mapped primary read를 HAP1/HAP2 양쪽 personal FASTA에 보낸다. 원래 unmapped read도 양쪽에 XU:i:1로 보존한다. HAP별 depth는 biological haplotype depth가 아니라 per-haplotype realignment evidence depth다.
- hp_routed_with_unknown: HP1은 HAP1, HP2는 HAP2로 보내고 HP0/missing 및 원래 unmapped read는 양쪽에 보존한다. coverage 보존용이며 HAP-specific label 용도가 아니다.
- hap_specific: 기존 HP tag의 일관성만 사용한다. AS 비교로 HAP을 새로 추정하지 않는다. Long read는 HP1/HP2 primary read만 각 HAP에 남긴다. Illumina는 QNAME fragment 단위로 HP를 모아 HP1-only fragment는 HAP1, HP2-only fragment는 HAP2에만 남긴다. HP0/missing, HP1/HP2 conflict, invalid HP fragment는 제외한다. 한 mate만 HP-tagged이면 그 mate의 HAP으로 pair를 보존한다. 원래 unmapped mate는 HP 증거 mate가 있는 fragment일 때만 같은 HAP에 XU:i:1로 보존한다.

각 sample QC 경로는 qc/<sample>.hap_read_routing.json이다. mode, 해석된 low-level policy, source/retained fragment 수, unknown/conflict/invalid-HP, one-mate-HP, unmapped-without-HP, HAP1/HAP2 QNAME 교집합을 기록한다. hap_specific은 하나의 mutually-exclusive dispatch branch를 사용하므로 QNAME 교집합은 0이어야 하며 0이 아니면 Stage 22가 실패해야 한다. Stage 22 summary TSV에도 hap_read_mode와 routing_qc_json을 기록한다.

기존 산출물은 mode-aware metadata 검사를 통과할 때만 재사용한다. 기존 all_to_both / hp_routed_with_unknown 메타데이터는 대응되는 mode로 호환한다. 과거 record-level strict 출력은 새 Illumina fragment-level hap_specific 출력으로 재사용하지 않는다.

## 260927 인수인계: 단일 진입점과 검증 완료 조건

- 단일 샘플: `PIPELINE_MANIFEST=... ./run_pipeline.sh sample <sample>`
- 목록: `PIPELINE_MANIFEST=... ./run_pipeline.sh sample-list <list> [count]`
- WGS 성공 뒤 manifest에 RNA 입력이 있는 샘플만 RNA가 실행된다. 샘플 하나의 실패는 목록의 다음 샘플을 막지 않는다.
- mode별 최종 경로는 `realigned/<mode>`와 `compact/<layout>__<mode>`이다. mode가 다른 결과를 같은 경로에 쓰지 않는다.
- Stage 12 기본은 `full` 내부 검증 후 독립 source trace까지 완료하는 것이다. WGS는 Stage 23, RNA pair는 Stage 25의 final package 검증을 통과한 뒤에만 각각 Stage 24/26 cleanup을 실행한다.
- RNA는 `qc/rna_preflight.json` fingerprint가 같은 산출물만 resume한다. RNA 상태는 `<sample>.rna.state.tsv`에 남는다.
- 카탈로그의 모델 사용 조건은 `training_eligible=True`다. `processed_qc_fail`은 기술 실패가 아니라 ASE QC 탈락이다.
- Clair3는 제공 VCF가 없을 때만 실행한다.
