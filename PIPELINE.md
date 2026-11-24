# original — WGS + RNA 전처리 단계별 설명

기준일: **2026-11-23**. GRCh38 + phased VCF로 HAP1/HAP2 서열을 생성하고 선택한 WGS/RNA payload를 저장함. 실행·JSON·자원·정리 설정은 [GUIDE.md](GUIDE.md)를 기준으로 함.

| WGS payload | Stage10 | Stage12 |
|---|---|---|
| `full` | 개인 HAP 재정렬 적용함 | HAP 서열·readset·HAP 기준 depth 저장함 |
| `sequence_only` | WGS-only/coverage_only에서 비적용임 | HAP 서열·window 주석·좌표 매핑 저장함 |
| `pileup_feature` | WGS-only/coverage_only에서 **비적용임** | reference 정렬 read pileup을 chain으로 HAP 좌표에 투영해 저장함 |

RNA `full`의 ASE/CNV용 WGS 근거는 WGS 출력 형식과 별개로 필요함. `MFM_PILEUP_BENCHMARK_REALIGN=1`은 실험용 Stage10 강제 옵션으로, 기본 pileup 처리시간에 포함하지 않음.
RNA13–26은 GTEx/TOPMed v11 기반 CPU 프로필(STAR 2.7.11b)을 사용함. 개인 HAP reference·lifted GENCODE v47·RNA window signal과 Stage20–23 ASE/CNV 근거는 이 파이프라인의 확장임. 실제 인자는 `rna_profile.py`와 실행 로그를 확인함.

## 실행

| 명령 | 하는 일 |
|---|---|
| `bash ./00_environment.sh` / `verify-core` / `verify-rna` | 인수 없으면 환경 설치·보강함. `verify-*`는 해당 도구 확인용임 |
| `./run_sample.sh <S>` | manifest 에 RNA 가 있으면 WGS → RNA, 없으면 WGS 만 |
| `./run_wgs_sample.sh <S>` | WGS 01–12 뒤 Stage 23 검증과 선택적 Stage 24 cleanup (`WGS_STOP_AFTER=NN`) |
| `./run_rna_sample.sh <S>` | RNA 13–26 (WGS 가 먼저 끝나 있어야 함) |
| `bash run_pipeline.sh dry-run` / `run --config FILE --executor local` | JSON·선택 TSV로 다운로드/계산/검증/정리를 이어가는 공식 내장 큐임 |
| `runner_local.sh`, `run_wgs_phased.sh`, `slurm_wgs.sbatch`, `slurm_rna.sbatch` | 여러 sample 을 로컬 병렬 또는 slurm array 로 |

설정: `config_v2.sh` (경로·도구·기본값), 서버 이식용 예시 `portable_config.example.sh`. sample 정보는 `manifests/` 에 있다: `preprocess_manifest.tsv`(3,212 sample 의 입력 정보), `read_files.tsv`(read 파일 목록), `portable_summary.json`(요약). 다른 manifest 는 `--manifest` 또는 `PIPELINE_MANIFEST` 로 지정한다.
`references/` 는 reference 파일로 가는 링크다. `references/rna_index_nonautosome_gtex.fa` 는 첫 RNA 실행 때 자동으로 만들어진다.

약칭: `F` = 최종 sample 폴더 `<output>/<S>` · `W` = sample work 폴더(WGS 24/RNA 26에서 삭제) · `R` = `F/rna_personal/<label>` · `h` = 1, 2

## WGS (01–12, Stage 23 검증, Stage 24 cleanup)

| # | 파일 | 하는 일 | 주 출력 |
|---|---|---|---|
| 01 | `01_prepare_analysis_bam.sh` | long-read uBAM/FASTQ 를 analysis reference 좌표 BAM 으로 표준화 (short read 는 BWA-MEM2, long read 는 minimap2) | 분석용 BAM |
| 02 | `02_preflight.sh` | sample 별 입력 계약 검사 (파일·index·reference·SM 등) | preflight QC |
| 03 | `03_prepare_vcf.sh` | sample VCF 준비. short read 는 1KGP phased panel 에서 추출, long read 는 제공 public VCF 검증·정규화 | seq / tag VCF |
| 04 | `04_call_longread_clair3.sh` | Clair3 self-call + WhatsHap. **기본 꺼짐** (`vcf_source=clair3` + `ENABLE_CLAIR3_SELF_CALL=1` 일 때만) | Clair3 VCF, HP/PS BAM |
| 05 | `05_clair3_qc.sh` | 04 산출물 QC | `qc/clair3_output_qc.*` |
| 06 | `06_haplotag.sh` | phased VCF 로 WhatsHap haplotag | HP/PS-tagged BAM |
| 07 | `07_denovo_assembly.sh` | (선택) hifiasm de novo assembly | assembly FASTA |
| 08 | `08_anchor_assembly.sh` | (선택) assembly contig 를 reference 에 anchoring QC. HAP FASTA 로 승격하지 않음 | anchoring QC |
| 09 | `09_personal_fasta.sh` | seq VCF(SNV + indel)로 HAP1/HAP2 개인 FASTA + reference↔HAP chain | `F/individual_fasta/HAP{1,2}/<chr>.fa·.fai·.chain`, `personal_fasta_qc.json` (모델 전용 정리 후에는 서열 Zarr·압축 chain·자산 manifest로 보존함) |
| 10 | `10_realign_wgs.sh` | full/RNA full 근거/실험 옵션에 적용함. WGS read를 HAP1/HAP2 개인 FASTA에 각각 재정렬. Illumina는 BWA-MEM2, long-read는 minimap2이며 HAP-specific은 명시 opt-in. mapped/paired fraction은 기본 record-only QC | `F/realigned/<hap_read_mode>/HAP{1,2}/<S>.personal.HAP{h}.cram` + mapped-fraction QC |
| 11 | `11_windows.sh` | **Stage10이 필요한 모드에서는 alignment 계약 확인 뒤**, 그 외에는 HAP 서열·chain 계약으로 1,048,576bp window index 생성. 비율 참고선은 strict일 때만 gate | `F/windows/windows.tsv` |
| 12 | `12_compact_zarr.sh` | full=readset/depth, sequence_only=서열, pileup_feature=서열+12채널 reference pileup임. 각 payload의 전수검증/독립 source trace를 수행함 | `compact/<payload-layout>/`, `qc/stage12_trace/*` |
| 23 | `23_validate_sample_package.sh` | WGS-only 또는 WGS–RNA full package의 sample ID/window/QC/log/provenance 완전성 검증 | `qc/final_package_validation.json`, `metadata/completeness.tsv` |
| 24 | `24_cleanup_intermediates.sh` | Stage 23 PASS일 때만 보수적 cleanup을 실행 | `metadata/cleanup_intermediates.tsv` |

## RNA (13–26)

```
RNA BAM ─(13)→ FASTQ ────────────────────────┐
GTF + chain ─(14)→ HAP GTF, exon interval ───┼─(15)→ STAR index ─(16)→ STAR BAM ─(17)→ dup 표시 BAM ─(18)→ coverage.zarr ─(19)→ QC JSON
                                             └── exon interval ──────────────┐
VCF + Stage 06 WGS BAM ─(20)→ WGS로 검증한 HAP variant ─────────────────────┐
                              HAP1/2 RNA BAM ────────────────────────────────┴─(21)→ allelic counts ─┐
                              WGS HAP CRAM + exon ─(22)→ gene rd ────────────────────────────────────┴─(23)→ ASE QC ─(24)→ sample.json ─(25)→ package 검증 ─(26)→ cleanup
```

아래 표는 RNA `full` 산출물의 의존 관계임. 실제 기본 실행은 HAP 병렬화와 RNA DAG를 사용할 수 있으며, 번호 순서가 곧 모든 작업의 직렬 순서를 뜻하지 않음. `coverage_only`는13–19와24–26을 실행하고20–23은 비적용으로 기록함. RNA 보류는 해당 branch를 실행하지 않음. `27_make_catalog.py`는 모든 sample 단계가 끝난 뒤 실행하는 독립 cohort utility다.

| # | 파일 | 하는 일 | 입력 | 출력 |
|---|---|---|---|---|
| 13 | `13_rna_to_fastq.py` | Stage 02에서 검증·고정한 RNA contract를 재확인한 뒤 BAM/CRAM → paired FASTQ | manifest RNA BAM/CRAM, Stage 02 contract | `W/rna/R1.fq`, `R2.fq`, `fastq_qc.json` |
| 14 | `14_lift_annotation.py` | GENCODE GTF 를 chain 으로 HAP 좌표로. 결실된 염기만 잘라내고, feature 전체가 결실일 때만 뺌 | `references/gencode.v47…gtf.gz`, HAP chain·`.fai` | `R/HAPh/annotation.gtf`, `exons.bygene.tsv`, `gtf_lift_qc.json` |
| 15 | `15_star_index.py` | HAP chr1–22 + GTEx index 의 나머지 contig(chrX·chrY·chrM·chrUn·random·EBV, alt·HLA·decoy 제외) + HAP GTF 로 STAR index. sjdbOverhang = 최대 read 길이 − 1 | HAP FASTA, HAP GTF, `fastq_qc.json` | `W/rna/idx_hap{h}_oh{N}_gtex_<GENCODE tag>_<fingerprint>/`, `hap{h}.index.json` |
| 16 | `16_star_align.py` | STAR 정렬 — GTEx `run_STAR.py` 설정 그대로 | FASTQ, index | `W/rna/aln_hap…/Aligned.sortedByCoord.out.bam`, `R/HAPh/Log.final.out`, `ReadsPerGene.out.tab` |
| 17 | `17_mark_duplicates.py` | Picard MarkDuplicates — GTEx `run_MarkDuplicates.py` 옵션, duplicate 는 표시만 → 최종 HAP BAM | STAR BAM | `R/HAPh/<S>.rna.personal.HAPh.bam(.bai)`, metrics |
| 18 | `18_coverage_zarr.py` | duplicate 포함 전장 1-bp coverage를 한 번 계산하고 WGS kept-window subset도 같은 실행에서 파생 | 17 BAM, HAP FASTA FAI, `F/windows/windows.tsv` | `R/HAPh/coverage.zarr`, `coverage_windows.zarr` |
| 19 | `19_rna_realign_qc.py` | HAP 별 정렬 QC, sample 요약. 재사용 판정(아래) | 13–18 JSON, coverage | `R/HAPh/rna_realign.json`, canonical `R/rna_realign_summary.json` (호환용 `rna_qc.json`) |
| 20 | `20_hap_variants.py` | HAP FASTA 생성 VCF의 exon ±1 kb 이형접합 SNV·indel을 Stage 06 WGS read로 확인. depth≥10, alt fraction≥0.1이면 `validated` | `personal_fasta_qc.json`의 VCF, haplotagged WGS BAM, reference, HAP FASTA/chain, GTF | `R/hap_variants.tsv`, `R/hap_variants_qc.json` |
| 21 | `21_allelic_counts.py` | 같은 fragment의 HAP1·HAP2 AS를 비교하되, 승리 HAP 정렬에서 Stage 20 validated variant를 덮는 경우만 계수. 미검증 variant를 덮으면 제외. duplicate 합집합 제외 유지 | HAP1/2 BAM, exon, `hap_variants.tsv` | `R/allelic_counts.tsv`와 sidecar, schema v5 meta |
| 22 | `22_cnv_filter.py` | WGS 개인 정렬 exon depth / 중앙값=rd. 두 HAP 모두 0.5–1.5이고 결측이 없어야 통과 | exon, Stage 10 CRAM | `R/gene_rd.tsv`와 policy meta |
| 23 | `23_ase_qc.py` | ASE 학습 적격성 및 v5 provenance 확인. gene별 validated/unvalidated variant 수와 제외 fragment 수 기록 | Stage 20–22 출력 | `R/ase_qc.tsv`, `ase_training_eligibility.tsv`, `ase_qc.json` |
| 24 | `24_stamp_rna_metadata.py` | RNA와 HAP variant/ASE 산출물 경로를 `sample.json`에 기록 | `F/` | `metadata/sample.json` |
| 25 | `25_validate_sample_package.sh` | WGS–RNA package와 Stage 20→21→23 schema/hash 연결, sample/window/QC/log/provenance 완전성 검증 | `F/` | `qc/final_package_validation.json`, `metadata/completeness.tsv` |
| 26 | `26_cleanup_intermediates.sh` | Stage25 PASS·현재 RNA payload의 근거가 모두 맞을 때만 정리함. coverage_only에 ASE 파일을 요구하지 않음 | `F/` | 모델 전용 설정에서 중간 BAM/CRAM/FASTA/index/work 정리함. 모델 입력·주석·압축 매핑·QC·metadata·log 보존함 |
| 27 | `27_make_catalog.py` | PASS sample cohort catalog. sample Stage 26 이후 batch 단위로 실행 | `F/metadata/sample.json` 들 | catalog TSV/JSON |

### RNA 설정과 재사용

- 정렬·coverage 설정(STAR·MarkDuplicates·coverage 인자, 추가 contig 규칙, lift 방식)은 `rna_profile.py` 한 곳에 있다.
- 이 설정과 GENCODE 파일 내용, 추가 contig 원본 `.fai` 내용으로 fingerprint 를 만든다. 다시 돌릴 때는 fingerprint 가 같을 때만 이전 결과를 쓴다.

| 단계 | 건너뛰는 조건 | 다시 만드는 경우 |
|---|---|---|
| 13 | `R1.fq`·`R2.fq`·`fastq_qc.json`이 있고 Stage 02 source fingerprint가 동일 | `RNA_FORCE=1` 또는 source contract 변경 |
| 14 | `gtf_lift_qc.json` schema 와 GTF 내용 sha256 이 같음 | GENCODE 파일 변경 |
| 15·16 | index·정렬 폴더가 있음 (폴더 이름에 sjdbOverhang·GENCODE tag·fingerprint) | read 길이·설정 변경 |
| 17 | marked BAM·BAI·metrics 가 있음 | `RNA_FORCE=1` |
| HAP 전체 (14–19) | 최종 BAM, `rna_realign.json`, 전장/Window coverage가 있으며 GENCODE tag·profile·Stage 02 contract가 같음 | 하나라도 다르면 |


## 보조 파일

| 파일 | 역할 |
|---|---|
| `src/integrate_wgs_rna/` | WGS 단계들이 부르는 Python 모듈 (VCF, 개인 FASTA, 재정렬, compact zarr 등) |
| `scripts/run_integrated_preprocess.sh`, `run_one_local.sh`, `resource_runner.py` | WGS 단계 실행기, 로컬 1-sample 실행, 자원 측정 |
| `rna_common.py` | RNA 13·15–17·19 공통 (로그, 외부 도구 실행, 원자적 JSON 쓰기) |
| `tools/personal_assets.py` | 영구 자산(HAP FASTA·chain) manifest 작성·검증 |
| `legacy_phase_longread.sh`, `anchor_blocks.py` | 예전 long-read phasing 경로 (기본 실행에서는 쓰지 않음) |
| `environment.yml` | conda 환경 |
