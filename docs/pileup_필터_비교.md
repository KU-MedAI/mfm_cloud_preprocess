# pileup과 readset의 전처리 정책 비교

확인일: 2026-11-24. 비교 대상은 현재 동봉된 `cloud_전처리용` 기반 코드의 `full` readset 모드와 `pileup_feature` 모드임. 모든 과거 실행의 저장 설정이 같았다고 단정하는 비교는 아님.

## 1. 판정

- 공통: HAP 서열 생성, window 선택, read의 MAPQ·SAM flag 제외 기준을 유지함.
- 변경: WGS readset 대신 reference에서 집계한 12채널 pileup을 HAP window에 투영해 저장함.
- 해석: **필터 기준이 같아도 정렬 대상이 다르므로 최종 read 집합과 depth가 같다고 보장할 수 없음.** 단순 파일 형식 변환으로 보지 않음.
- 현재 생물학 소스: `src/integrate_wgs_rna`와 `vendor`의 배포 파일175개가 점검 시점의 `original`과 바이트 단위로 동일함. 이번 준비에서는 필터·정렬·feature 수식을 수정하지 않음.

## 2. read 필터의 실제 적용값

두 모드는 [공통 프리셋](../src/integrate_wgs_rna/full_preprocess.py)의 `resolve_stage70_filter_policy()`를 호출함. Python의 record 필터는 [rejection_reason()](../vendor/MFM-DNARNAFM/pre_training/preprocess_raw_readset/filters.py)를 공유함. pileup의 C++ 생산기/독립 검증기도 같은 설정에서 제외 flag와 MAPQ를 전달받음.

| 항목 | 현재 readset full | 현재 pileup_feature |
|---|---|---|
| MAPQ | 기본30 이상 | 동일. `MIN_MAPQ`를 전달함 |
| Unmapped | 제외 | 동일 |
| Secondary | 제외 | 동일 |
| QC fail | 제외 | 동일 |
| read sequence 없음 | 제외 | 동일 |
| Short read duplicate | 제외 | 동일 |
| Short read supplementary | 제외 | 동일 |
| Long read duplicate | 보존 | 동일 |
| Long read supplementary | 보존하되 주 모델·depth는 primary만 사용하고 supplementary는 보조 track으로 분리함 | 주 feature·depth는 primary만 사용하고 `supplementary_depth`로 분리함 |
| Proper pair | 기본 강제하지 않음 | 동일 |
| NM·mismatch·길이·clipping 비율 | 기본 추가 read 제외 기준 없음 | 동일 |
| Primary base depth | CIGAR `M/=/X`를 집계함. `D/N/I/S/H/P`를 base depth에 더하지 않음 | 동일한 base 관측 규칙임. 다만 좌표는 reference 기준임 |

**현재 Long read 기본 정책은 duplicate 전체 제외가 아님.** `long_read_primary_supp_evidence_v1`이 적용되며, Short read의 `short_read_base_v1`과 구분됨. 이는 pileup으로 변경하면서 새로 바꾼 규칙이 아님.

- BQ: pileup은 낮은 BQ를 이유로 염기 관측을 제거하지 않음(`min_base_quality=0`). BQ가 없거나255이면 depth에는 남기고 평균 BQ의 분모 `baseq_depth`에서는 제외함.
- BAQ: pileup에서 재계산하지 않음. mate overlap: 각 정렬을 따로 집계함. depth cap: 없음. samtools mpileup의 기본 BQ/중복 mate/max-depth 옵션을 그대로 적용하는 구현이 아님.
- 서열의 N: read 안의 N을 이유로 read 전체를 제외하지 않음. primary 관측의 `read_N_ratio`로 집계함. 지원하지 않는 관측 염기 문자는 오류로 처리함.
- 검증: [필터 설정 동등성 검사](../tests/test_pileup_feature.py), [native/Python·독립 HTSlib 검사](../tests/test_pileup_native.py)에서 MAPQ 경계·flag·BQ·CIGAR·9,001 depth 등을 확인함.

## 3. 실제 계산에서 달라지는 부분

| 항목 | readset full | pileup_feature |
|---|---|---|
| WGS 관측 정렬 | 개인 HAP1·HAP2 재정렬 결과 | Stage06 reference haplotag BAM |
| 저장 | HAP별 read sequence·정렬/CIGAR 메타데이터와 depth | HAP 서열과 reference 기반 12채널 feature·보조 배열 |
| HAP1/2 관측 구분 | 기본 `all_to_both`로 각각 재정렬한 결과를 사용함 | 동일 reference 위치의 관측을 두 HAP에 투영함. HP 태그로 나눈 allele별 pileup이 아님 |
| HAP 삽입 위치 | HAP 재정렬 read의 실제 염기 관측을 사용 가능함 | 삽입 DNA는 그대로 보존함. reference 대응이 없어 base ratio·depth·BQ는 만들지 않고 왼쪽 anchor의 generic insertion 비율만 전달함 |
| HAP 삭제 위치 | 개인 HAP의 존재하는 위치에서 관측함 | HAP에 없는 reference 위치는 HAP 입력의 염기 위치로 저장되지 않음 |
| read 복원 | 저장 read sequence를 읽을 수 있음 | 요약 tensor에서 개별 read·read ID·원래 염기 배열을 복원할 수 없음 |

reference에서 MAPQ29였던 read가 HAP 재정렬에서 MAPQ60이 되거나 그 반대가 되는 경우, 두 모드의 MAPQ30 기준은 같아도 선택 결과가 달라질 수 있음. primary/supplementary 판정·CIGAR도 정렬에 따라 바뀔 수 있음.

pileup 삽입 위치의 `depth=0`은 관측 read가 실제로 전혀 없다는 뜻과 구분해야 함. `mapping_kind=2`는 reference 비대응 삽입임. `mapping_kind=0`은 padding,1은 reference 대응 위치임.

12채널은 `log1p_depth`, `mean_mapq`, `mean_baseq`, `strand_forward_ratio`, `read_A/C/G/T/N_ratio`, `insertion_ratio`, `deletion_after_event_ratio`, `deletion_after_max_len`임. 삽입 allele 일치율·삽입 각 염기의 관측은 산출하지 않으며 DNA 서열과 삽입 사건 지지를 함께 제공함.

## 4. window·N·blacklist는 같은가

두 모드 모두 [11_windows.sh](../11_windows.sh) → [window_index.py](../src/integrate_wgs_rna/window_index.py)의 같은 window 원장을 사용함. pileup export가 새 window를 재선별하지 않음.

| 항목 | 제공 기본 설정에서의 공통 정책 |
|---|---|
| 염색체 | chr1–22 |
| HAP 서열 | GRCh38 + phased SNV/indel VCF의 consensus |
| Window | 1,048,576bp·비중첩·seed42의 sample/HAP/염색체별 결정적 시작점 |
| 말단 N | 각 염색체/HAP 양끝 연속 N만 trim |
| 내부 N | 보존. 실제 구간에 A/C/G/T가 하나도 없는 창만 제외 |
| Blacklist | 겹침 위치·비율 기록. 기본 `WINDOW_MAX_BLACKLIST_FRAC=1.0`으로 제거하지 않음 |
| N 기록 | `n_count`, `n_nonacgt`, `n_fraction`, `n_max_run` 등 원값을 기록함. 1%/5%를 저장 필터로 적용하지 않음 |
| 짧은 마지막 창 | 실제 길이를 기록하고 나머지는 P로 padding함 |
| 좌표 대응률 | 기록만 하며 window 제외 기준으로 사용하지 않음 |

`n_fraction`의 실제 정의는 `n_nonacgt / actual_length`임. N을 포함한 비 ACGT 염기 비율이며 padding은 분모에서 제외함. A/C/G/T/N 서열에서는 일반적인 N 비율과 같음.

사용자가 `WINDOW_DROP_N=1` 또는 더 낮은 `WINDOW_MAX_BLACKLIST_FRAC`를 명시하면 선택이 달라질 수 있음. 현재 제공 pileup 예제는 그러한 과거의 엄격 필터를 켜지 않음. 작은 통합검사에서 pileup 생성 전후 원본 `windows.tsv`의 바이트가 같음을 확인함.

## 5. 개인 HAP 재정렬·RNA 처리 조건

[config_v2.sh](../config_v2.sh)의 `sample_hap_realign_required()`와 [run_wgs_sample.sh](../run_wgs_sample.sh)가 결정함.

| 출력 설정 | WGS 개인 HAP 재정렬 | 이유 |
|---|---|---|
| readset full | 수행 | HAP별 readset 생성 |
| pileup + WGS-only | 기본 비적용 | reference 관측과 consensus chain으로 feature 생성 |
| pileup + RNA coverage_only | 기본 비적용 | RNA expression까지만 생성하며 ASE/CNV는 비적용 |
| pileup + RNA full | 수행 | 기존 RNA ASE/CNV가 요구하는 개인 WGS 근거 준비 |
| pileup + 비교 실험 옵션1 | 수행 | `MFM_PILEUP_BENCHMARK_REALIGN=1`일 때의 별도 실험임. pileup 관측은 여전히 reference BAM 사용 |

- 현재 Docker·cloud pileup 예제: `MFM_RNA_OUTPUT_MODE=full`, `MFM_SEQUENCE_DEFER_RNA=0`, `MFM_PILEUP_BENCHMARK_REALIGN=0`임. 따라서 pair의 ASE/CNV도 포함하며 WGS HAP 재정렬이 실행될 수 있음.
- RNA 정렬·1bp expression·ASE/CNV의 계산 및 read 기준은 pileup 출력 모드 때문에 바뀌지 않음. `coverage_only`를 별도로 선택하면 분석 범위가 달라지므로 full과 구분해야 함.
- 모든 모드에서 최종 QC와 source trace 후 정리함. pileup은 reference BAM에서 독립 재집계해 tensor를 대조하며, readset은 HAP 정렬/readset 계약을 검증함. 입력 좌표에 맞춰 검증 대상이 달라짐.

## 6. 코드만으로 실행 가능한 범위

- 포함: 전처리 코드와 vendor 모듈, 자동 환경 설치, WGS/RNA 단계, 입력 목록·공통 자산 준비, 내장 큐, QC·재개·정리임. 빌드할 때 다른 전처리 저장소의 소스 접근은 필요하지 않음.
- 실행 조건: Linux x86_64, Docker 또는 기본 OS 도구, 입력 URL 네트워크, 쓰기·잠금·rename·native 실행을 지원하는 디스크, CPU/RAM 설정이 필요함. 원본/reference 실물은 이미지에 넣지 않음.
- CPU Docker 범위: 제공 public/panel VCF 기반 WGS/RNA를 대상으로 함. GPU Parabricks·CuPy·Clair3 self-call은 현재 CPU 이미지에 포함하지 않음. 일반 설치의 Clair3는 별도 환경 설치 경로가 있음.
- 검증 범위: 작성 서버의 환경·소형검사와 GitHub 자동 검사의 실제 image build/run 결과를 구분함. 커밋별 PR Checks에서 이미지 검사 상태를 확인함. 새 VM의 대용량 donor 전과정 완주는 별도 검증 대상임.
