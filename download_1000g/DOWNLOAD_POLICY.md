# Download and preprocessing policy

> **현재 적용 범위:** 이 정책은 2026-08의 1KGP 3,202명 패키지 기록이다. current production cohort·VCF·checksum 정책은 `../share_package/DOWNLOAD_POLICY.md`가 우선한다.

## 목적

이 문서는 `samples_3202.csv`, `samples_3202.tsv`, `phase_assignments.tsv`를 이용해 1000 Genomes 3,202 sample 후보를 분할 다운로드하고 전처리하는 운영 정책임.
모든 경로는 이 package 폴더 기준 상대경로임.

## 핵심 정책

- 공용 resource는 phase 0으로 분리함.
- sample 데이터는 phase 1-5로 나눔.
- phase는 sample 수가 아니라 예상 raw 다운로드 용량 기준으로 균등화함.
- 각 phase는 다시 shard로 나눠 다운로드함.
- 다운로드 없이 전처리만 수행할 수 있어야 함.
- CSV는 사람 확인용, TSV는 bash/awk 자동화용으로 사용함.
- recount3 BigWig과 non-BAM RNA source는 이번 다운로드/전처리 운영 대상에서 제외함.

## 기본 실행 단위

| 단위 | 의미 | 파일 |
|---|---|---|
| phase 0 | common VCF/reference 등 공용 데이터 | `scripts/04_download_common_vcf.sh` |
| phase 1-5 | sample-level 논리 배치 | `phase_assignments.tsv` |
| shard | phase 내부 실제 다운로드 단위 | `download_lists/phase_XX/shards/` |

## 분할 정책

기본값은 다음과 같음.

```bash
PHASES=5
BALANCE_MODE=aligned_ready_fallback
INCLUDE_RNA=1
INCLUDE_BIGWIG=0
SHARD_GIB=1500
```

`BALANCE_MODE`는 다음 중 하나를 선택함.

| mode | 의미 | 권장 상황 |
|---|---|---|
| `aligned_ready_fallback` | Illumina/ONT Vienna aligned CRAM은 실제 크기 사용, HPRC/HiFi raw는 fallback estimate 사용 | 기본 운영 |
| `primary_url` | sample별 대표 URL 1개 기준 | 빠른 pilot 또는 URL smoke test |
| `all_manifest_runs` | manifest의 sample/source 전체 run 기준 | raw multi-run 전체를 실제로 받을 때 |

## source 선택 정책

`WGS_SOURCE_MODE`는 URL list 생성 시 사용함.

| mode | 의미 |
|---|---|
| `assigned` | `samples_3202.tsv`의 assigned source와 primary URL을 사용함. 기본값임. |
| `primary` | 대표 URL만 사용함. multi-run coverage가 부족할 수 있음. |
| `all_manifest_runs` | `manifests/wgs_urls.tsv`에서 sample/source에 해당하는 모든 run URL을 사용함. |
| `aligned_only` | Illumina 30x와 ONT Vienna aligned CRAM만 URL에 포함함. HPRC raw/HiFi raw는 제외함. |

`RNA_SOURCE_MODE`는 다음 중 하나임.

| mode | 의미 |
|---|---|
| `strict_bam` | GEUVADIS TOPMed RNA BAM만 포함함. 기본값임. |
| `primary` | `samples_3202.tsv`의 대표 RNA URL만 포함함. |
| `none` | RNA URL을 만들지 않음. |

## 추천 운영

### 1. 공용 VCF 준비

```bash
DRY_RUN=1 bash scripts/04_download_common_vcf.sh
bash scripts/04_download_common_vcf.sh
```

### 2. 기본 5 phase 생성

```bash
bash scripts/00_make_phase_assignments.sh
```

### 3. phase 1 URL list 생성

```bash
SHARD_GIB=1500 INCLUDE_RNA=1 INCLUDE_BIGWIG=0 bash scripts/01_make_phase_url_lists.sh 1
```

### 4. phase 1 shard 1 다운로드

```bash
DRY_RUN=1 bash scripts/02_download_phase.sh 1 shard_001
bash scripts/02_download_phase.sh 1 shard_001
```

### 5. 검증

```bash
bash scripts/03_verify_phase.sh 1 shard_001
```

### 6. 다운로드 없이 전처리만 수행

이미 raw가 준비된 경우 다운로드 단계 없이 sample list만 생성하고 전처리 command를 주입함.

```bash
MODE=pre PLATFORM=ONT bash scripts/05_preprocess_phase.sh 1 shard_001
```

실제 전처리 실행 예시는 다음과 같음.

```bash
RUN_PREPROCESS=1 \
PROCESS_CODE_DIR=../process_code \
PREPROCESS_CMD='SAMPLES="$SAMPLES_SPACE" ./run_cohort.sh igsr_ont_wgs_only' \
bash scripts/05_preprocess_phase.sh 1 shard_001
```

## 주의사항

- HiFi/CCS와 HPRC ONT는 multi-run일 수 있음. `primary` 모드는 coverage가 부족할 수 있음.
- HGSVC3 manifest에는 `Kinnex`/`*.flnc.bam` RNA 성격 파일이 섞일 수 있음. WGS URL 생성 스크립트는 `HGSVC_HiFi` source에서 해당 파일을 제외함.
- `all_manifest_runs`는 용량이 매우 커질 수 있으므로 pilot 없이 전체 실행하지 않음.
- strict RNA-pair는 RNA BAM 기준 449개임. RNA BAM이 없는 sample은 모두 WGS-only로 통합함. BigWig과 non-BAM RNA source는 다운로드 대상에서 제외함.
- 현재 package는 다운로드와 전처리 운영을 분리함. 다운로드 없이 바로 전처리하려면 `05_preprocess_phase.sh`를 사용함.
