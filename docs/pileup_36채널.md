# 기본 pileup 36채널과 해석

2026-11-25 기본 출력임. 과거 raw 기본30(공통12+HP1별9+HP2별9)에 공통6개를 추가함. 동일 지표를 중복하지 않으며99채널 구성은 사용하지 않음.

## 실행·출력

- Short/Long 모두 미지정 기본값은 `pileup_feature`임. 기존 readset은 `full`, 서열만은 `sequence_only`로 명시함. RNA는 기존 `full`이 기본임.
- 저장: `compact/hap_sequence_reference_pileup_v2/<chr>/pileup.zarr/HAP1|HAP2/features`임. 형상은 `[padded position, 36]`, float32·zstd level5임. HAP DNA/window·read 필터·RNA/ASE/CNV 정책은 유지함.
- schema는 `mfm_reference_pileup_projected_v2`, feature_set은 `reference_hp_baf_v2_36ch`임. 실제 채널 순서는 manifest와 아래 표에 기록함.
- 같은 reference 위치에 대응하는 두 HAP window는 같은 reference 관측을 받음. 각 window 안에 공통/HP1/HP2 채널이 모두 존재함. `hp1_`/`hp2_`는 입력 BAM의 정수 HP=1/2 read 집계이며 window HAP와 혼동하지 않음. HP가 없거나 유효하지 않은 read는 공통에만 기여함. HP와 개인 HAP명은 제공된 phase/VCF 근거를 따르며 새로 정렬·추정하지 않음.

## 분모·결측

- 기본 depth: primary CIGAR M/=/X에서 관측한 A/C/G/T/N 합계임. deletion span은 base depth에서 제외함. HP9도 해당 HP의 같은 분모를 사용함.
- `mean_mapq`: 공통 base depth 분모임. `mean_baseq`: 실제 BQ가 있는 baseq_depth 분모임. missing QUAL/255를 BQ0 관측으로 만들지 않음.
- forward/reverse: 공통 base depth에 대한 각 방향의 관측 비율임. depth>0일 때 두 비율의 합이1임. 공통 forward는 기존 채널을 유지하고 reverse만 추가함.
- insertion/deletion 사건 비율: 왼쪽 anchor의 base depth+deletion span depth가 분모임. 최대 deletion 길이는 관측 read CIGAR의 길이임.
- BAF3: 개인 sequence VCF의 유효 biallelic A/C/G/T SNV·diploid GT만 사용함. `vcf_ref_allele_fraction=REF/(REF+ALT)`, `baf=ALT/(REF+ALT)`, `minor_allele_fraction=min(두 비율)`임. N/다른 allele은 이 분모에 포함하지 않음. 다중 ALT·indel·missing GT·동일 위치 중복은 BAF site에서 제외하되 window/read는 제거하지 않음. VCF REF 불일치는 실패로 처리함.
- `read_supports_variant_ratio`: reference와 다른 A/C/G/T 중 가장 많은 염기 관측수/base depth임. `base_diff_from_ref_ratio`: N을 포함해 reference와 다른 관측 전체/base depth임. 둘 다 A/C/G/T reference에서 정의하며 VCF 지정 ALT의 BAF와 다름.
- `baf_site`: 지원하는 VCF site임. `baf_heterozygous`: 그 site의 GT가 이형접합임. `baf_valid`: site이며 REF+ALT 관측 분모가 양수임. 숫자0과 미정의 위치를 mask로 구분함. HP별 BAF를 별도로 만들지 않음.
- HAP 삽입 token: 기존 chain의 anchor에서 공통/HP1/HP2 insertion_ratio만 전달함. 삽입 각 염기의 depth·BQ·BAF를 생성하지 않음. reference 대응이 없는 삽입/padding의 reference_base_code는255임. reference 대응 위치는 A/C/G/T/N=0/1/2/3/4임.

## 채널 순서

| 0부터 순서 | 채널 | 분류 |
|---:|---|---|
| 0 | `log1p_depth` | 공통 기본임 |
| 1 | `mean_mapq` | 공통 기본임 |
| 2 | `mean_baseq` | 공통 기본임 |
| 3 | `strand_forward_ratio` | 공통 기본임 |
| 4 | `read_A_ratio` | 공통 기본임 |
| 5 | `read_C_ratio` | 공통 기본임 |
| 6 | `read_G_ratio` | 공통 기본임 |
| 7 | `read_T_ratio` | 공통 기본임 |
| 8 | `read_N_ratio` | 공통 기본임 |
| 9 | `insertion_ratio` | 공통 기본임 |
| 10 | `deletion_after_event_ratio` | 공통 기본임 |
| 11 | `deletion_after_max_len` | 공통 기본임 |
| 12 | `hp1_log1p_depth` | HP1 기본임 |
| 13 | `hp1_read_A_ratio` | HP1 기본임 |
| 14 | `hp1_read_C_ratio` | HP1 기본임 |
| 15 | `hp1_read_G_ratio` | HP1 기본임 |
| 16 | `hp1_read_T_ratio` | HP1 기본임 |
| 17 | `hp1_read_N_ratio` | HP1 기본임 |
| 18 | `hp1_insertion_ratio` | HP1 기본임 |
| 19 | `hp1_deletion_after_event_ratio` | HP1 기본임 |
| 20 | `hp1_deletion_after_max_len` | HP1 기본임 |
| 21 | `hp2_log1p_depth` | HP2 기본임 |
| 22 | `hp2_read_A_ratio` | HP2 기본임 |
| 23 | `hp2_read_C_ratio` | HP2 기본임 |
| 24 | `hp2_read_G_ratio` | HP2 기본임 |
| 25 | `hp2_read_T_ratio` | HP2 기본임 |
| 26 | `hp2_read_N_ratio` | HP2 기본임 |
| 27 | `hp2_insertion_ratio` | HP2 기본임 |
| 28 | `hp2_deletion_after_event_ratio` | HP2 기본임 |
| 29 | `hp2_deletion_after_max_len` | HP2 기본임 |
| 30 | `vcf_ref_allele_fraction` | 공통 BAF임 |
| 31 | `baf` | 공통 BAF임 |
| 32 | `minor_allele_fraction` | 공통 BAF임 |
| 33 | `read_supports_variant_ratio` | 공통 보완임 |
| 34 | `base_diff_from_ref_ratio` | 공통 보완임 |
| 35 | `strand_reverse_ratio` | 공통 보완임 |

## 보조 배열·검증·자원

- uint32: 공통/HP1/HP2의 depth·supplementary_depth·max_deletion_length와 공통 baseq_depth임. uint8: mapping_kind·reference_base_code·BAF mask3개임. offsets/lengths는 기존 window와 같음.
- 생산: native CIGAR 순회 한 번으로 공통과 해당 HP 관측을 동시에 계산함. VCF/reference는 block별 조회하며 HAP1/2가 캐시를 공유함. production block262,144bp·cache512MiB, source 검사 block131,072bp임. 전체 염색체 tensor를 RAM에 올리지 않음.
- worker 기본 RAM예약4GiB·주소공간 보호8GiB임. cache512MiB + 1M window tensor약144MiB + 보조 배열/원시 count·decode·chain·압축 버퍼를 포함하며 샘플 CPU/RAM 예산에 맞춰 병렬도를 정함. 생산 cache를 반환한 후 독립 검증을 수행함. 극단적인 입력은 필터를 완화하지 않고 RAM 예산을 재검토함.
- 검증: 독립 HTSlib pileup과 저장된 모든 배열을 전수 대조함. source/VCF/sample/policy 지문과 SHA·형상·mask를 확인함. v2 source 검증은 native HTSlib가 필요하며 빌드 불가 때 검사 생략으로 완료하지 않음. 생산만 compiler가 없을 때 동등한 Python 경로를 지원함.
- 구버전12채널은 reader·검사에서 지원함. 신규36채널 완료 판정은 실제 manifest와 최종QC의 feature_set이 일치해야 함. 구버전 결과를 새 완료로 건너뛰거나 덮어쓰지 않음. 새 처리는 별도 output/state/run_id로 실행함.

## 주요 코드

- 채널·저장 계약: `src/integrate_wgs_rna/pileup_features/schema.py`임.
- 단일 순회: `native/pileup_counts_complete.cpp`, `native/pileup_complete.hpp`임. 독립 집계: `native/pileup_oracle_complete.cpp`임.
- BAF site: `pileup_features/variants.py`임. 계산·bounded cache: `complete_compute.py`임. 저장·source 검증: `complete_export.py`임.
- 공식 진입: `11_windows.sh`·`12_compact_zarr.sh` → `pileup_features/export.py`임. prepared sequence VCF와 실제 source sample ID를 전달함.
- 읽기: `WindowReader.load()`의 `pileup_features`, `pileup_channels`, HP depth·BAF mask를 사용함. 기존 DNA·RNA 필드는 유지함.
