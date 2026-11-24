# 1000 Genomes 3,202 sample share package

> **상태(2026-09-16): historical 1KGP downloader.** 이 폴더는 2026-08의 3,202명·ONT 포함 배정 기록과 1KGP 다운로드 도구를 보존한다. 현재 PacBio 우선 3,212명 registry와 전달 문서는 `../share_package/`를 사용한다. 이 폴더의 공용 1KGP VCF 정책을 PacBio 295명에 적용하면 안 된다.

## 개요

이 폴더는 1000 Genomes 30x Illumina WGS가 존재하는 3,202명에 대해 WGS-only/RNA-pair/platform별 사용 방식을 정리한 전달용 패키지임. 모든 내부 경로는 이 폴더 기준 상대경로로 작성함. 서버 local 절대경로는 포함하지 않음.

## 파일 구성

| file | description |
|---|---|
| `samples_3202.csv` | sample-level 분류표와 source URL. |
| `summary_counts.csv` | 주요 집계표. |
| `summary.json` | 주요 집계의 JSON 요약. |
| `source_projects.md` | 원천 project, 공식 링크, 다운로드 예시. |
| `DOWNLOAD_GUIDE.md` | 데이터 형식별 다운로드 방법과 검증 명령. |
| `DOWNLOAD_POLICY.md` | 5 phase/shard 분할, source 선택, no-download 전처리 정책. |
| `scripts/` | phase 생성, URL list 생성, 다운로드, 검증, 전처리 wrapper bash 스크립트. |
| `manifests/` | sample/source별 다운로드 manifest와 분류 근거 파일. |

## 최종 추천 사용 그룹: strict RNA BAM 기준

| recommended_training_use | assigned_wgs_platform | sample_count |
|---|---|---:|
| post_training_wgs_rna_pair | HiFi_CCS | 28 |
| post_training_wgs_rna_pair | Illumina_30x | 304 |
| post_training_wgs_rna_pair | ONT | 117 |
| pre_training_wgs_only | HiFi_CCS | 215 |
| pre_training_wgs_only | Illumina_30x | 1700 |
| pre_training_wgs_only | ONT | 838 |

해석: RNA BAM이 확인된 449개만 post-training WGS-RNA pair 후보임. RNA BAM이 없으면 RNA raw/FASTQ 흔적 여부와 무관하게 모두 WGS-only pre-training으로 통합함.

## 전체 WGS platform priority 기준

| priority_group | sample_count |
|---|---|
| 01_HiFi_CCS_long_read | 243 |
| 02_ONT_long_read | 955 |
| 03_Illumina_short_read | 2004 |

우선순위는 `HiFi/CCS > ONT > Illumina`임. 같은 sample이 여러 WGS platform을 가지면 가장 높은 우선순위 한 그룹에만 배정함.

## strict RNA BAM priority 기준

| priority_group | sample_count |
|---|---|
| 01_RNA_BAM_plus_HiFi_CCS | 28 |
| 02_RNA_BAM_plus_ONT_no_HiFi | 117 |
| 03_RNA_BAM_plus_Illumina_remainder | 304 |
| 04_HiFi_CCS_long_read_no_RNA_BAM | 215 |
| 05_ONT_long_read_no_RNA_BAM | 838 |
| 06_Illumina_short_read_no_RNA_BAM_longread | 1700 |

RNA-pair는 RNA BAM이 있는 경우만 인정함. RNA BAM이 없으면 모두 WGS-only로 처리함.

## WGS source별 URL coverage

| wgs_download_source | sample_count |
|---|---|
| ENA_PRJEB31736_1kGP_30x | 2004 |
| HGSVC_HiFi | 17 |
| HPRC_HiFi | 196 |
| HPRC_ONT | 30 |
| ONT_Vienna | 863 |
| missing_or_lookup_required | 92 |

`missing_or_lookup_required`는 donor matrix에는 platform flag가 있지만 package에 포함된 URL manifest에서 해당 sample의 URL을 찾지 못한 경우임.


## 운영 스크립트

- `scripts/00_make_phase_assignments.sh`: 5 phase 또는 사용자 지정 phase 수로 sample을 균등 분할함.
- `scripts/01_make_phase_url_lists.sh`: phase/shard별 다운로드 URL list를 생성함.
- `scripts/02_download_phase.sh`: phase/shard URL list를 `aria2c` 또는 `wget`으로 다운로드함.
- `scripts/03_verify_phase.sh`: 다운로드된 CRAM/BAM과 index를 검증함.
- `scripts/04_download_common_vcf.sh`: 2022 1kGP phased VCF chr1-22를 다운로드함.
- `scripts/05_preprocess_phase.sh`: 다운로드 없이 phase/shard sample list를 만들고, 필요 시 사용자가 지정한 전처리 명령을 실행함.

상세 정책은 `DOWNLOAD_POLICY.md`, 실제 명령 예시는 `DOWNLOAD_GUIDE.md`에 있음.

## 다운로드 원칙

- 실제 다운로드 방법과 데이터 형식은 `DOWNLOAD_GUIDE.md`에 정리함.
- sample별 다운로드는 `samples_3202.csv`의 `wgs_primary_url`, `rna_primary_url`을 먼저 확인함.
- multi-run source는 primary URL 하나만 받으면 부족할 수 있음. 반드시 `wgs_source_manifest` 또는 `rna_source_manifest`에서 같은 sample/source의 전체 URL을 필터링함.
- 공용 VCF는 모든 sample에서 같은 2022 1kGP high coverage phased panel을 사용함.
- 원천 project와 공식 root link는 `source_projects.md`에 정리함.

## 주요 CSV 컬럼

- `recommended_training_use`: 현재 추천 용도. `post_training_wgs_rna_pair` 또는 WGS-only 계열임.
- `wgs_pair_type_strict`: `RNA-pair` 또는 `WGS-only` 간단 구분임.
- `assigned_wgs_platform`: 중복 제거 후 사용할 WGS platform임.
- `assigned_wgs_project`: 해당 WGS 데이터의 원천 project명임.
- `assigned_wgs_project_root_url`: 원천 project/root 링크임.
- `wgs_primary_url`: 대표 WGS 다운로드 URL임.
- `wgs_source_manifest`: 전체 WGS URL을 확인할 package-relative manifest임.
- `rna_project_available`: RNA 관련 데이터의 원천 project명임.
- `rna_primary_url`: 대표 RNA 다운로드 URL임.
- `rna_source_manifest`: 전체 RNA URL을 확인할 package-relative manifest임.
- `rna_bigwig_url_validated_subset`: 과거 audit 흔적으로 남긴 recount3 BigWig URL임. 이번 다운로드/전처리 대상에서는 제외함.
- `common_vcf_panel_url_template`: 2022 1kGP phased panel VCF URL template임.

## 주의사항

- 이 패키지는 sample 분류와 다운로드 출처 전달용임. 실제 대량 다운로드 전에는 각 URL의 존재 여부와 checksum을 다시 확인하는 것이 안전함.
- HPRC HiFi/ONT는 raw uBAM/BAM run일 수 있어 현재 전처리 파이프라인에 바로 넣기 전에 alignment 여부를 확인해야 함.
- RNA BAM이 없는 sample은 RNA 후보 여부를 따로 분류하지 않고 WGS-only로 통합함. recount3 BigWig은 이번 다운로드 대상에서 제외함.
