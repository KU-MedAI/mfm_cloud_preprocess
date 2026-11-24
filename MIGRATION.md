# v2 적용 메모 (2026-09-25)

> 아래 표는 과거 번호를 바꾼 migration 기록이며, 그 표의 RNA 번호가 현재 번호라는 뜻은 아니다. v2 현재 실행 번호는 WGS `01–12`(선택 assembly 분기 `07–08` 포함), RNA `13–26`, cohort catalog `27`이다. 실제 호출은 `run_wgs_sample.sh`, `run_rna_sample.sh`, `run_pipeline.sh`를 기준으로 확인한다.

> 중요: 구번호 파일 중 일부는 얇은 wrapper가 아니라 해당 시기의 계산을 직접 수행하는 독립 구현이다. 파일을 보존했으며, 일부 helper만 `rna_*_common.py`로 공유한다. 구번호 진입점을 실행하면 현행 단계와 다른 입력 요건이나 통계 정책을 쓸 수 있으므로 새 전처리에는 사용하지 않는다.

## 현재 v2 실행 단계와 구번호 파일

| 구번호 파일 | 현재 v2 단계 | 상태와 차이 |
|---|---|---|
| `13_lift_annotation.py` | Stage 14 | 구 실행 파일. 현재 Stage 13은 RNA FASTQ 복원, Stage 14는 annotation lift다. |
| `14_realign_rna_star.sh` | Stage 16–17 | 구 통합 STAR 단계. 현재는 Stage 16 정렬과 Stage 17 MarkDuplicates로 나뉜다. |
| `15_coverage_zarr.py` | Stage 18 | 구 coverage 구현. 현재 runner는 `18_coverage_zarr.py`를 사용한다. |
| `16_allelic_counts.py`, `20_allelic_counts.py` | Stage 21 | 구 allele count 구현. 현행 Stage 21은 Stage 20의 검증 variant 증거를 요구한다. |
| `17_cnv_filter.py`, `21_cnv_filter.py` | Stage 22 | 구 CNV 진입점. 현행 Stage 22는 현재 RNA 단계의 병렬·QC 계약을 따른다. |
| `18_ase_qc.py`, `22_ase_qc.py` | Stage 23 | 구 ASE QC. 현행 Stage 23은 검증된 HAP variant 기반 입력과 별도 schema를 요구한다. |
| `19_stamp_rna_metadata.py`, `23_stamp_rna_metadata.py` | Stage 24 | 구 metadata 진입점. 공통 stamp 동작만 공유하며 Stage 24는 더 많은 QC 경로를 기록한다. |
| `21_make_catalog.py` | Stage 27 | 이전 catalog schema용 독립 구현. 현재 catalog 생성기는 `27_make_catalog.py`다. |

v2에는 개인 assembly 분기 Stage 07–08이 있어 RNA 번호가 `to_tgil` 전달본과 같지 않다. 위 현재 단계 번호는 v2 runner를 따르며 다른 파이프라인의 번호로 맞추지 않는다.

WGS 구현은 BWA-MEM2·chain·legacy Stage70 검증 schema를 보존한다. Python schema·CLI의 `stage22`/`stage70` 명칭은 호환성을 위해 유지하지만 현재 실행 단계는 각각 Stage 10/12다.

# 스크립트 재번호 기록 (2026-09-17; 당시 snapshot)

번호를 **실행 순서**에 맞췄다. 이전 번호는 순서와 무관해서 읽기 어려웠다:

- `22_realign_wgs.sh`가 `41_windows.sh` **뒤에** 돌았다
- `20_`으로 시작하는 파일이 두 개였다 (`20_call_longread_clair3.sh`, `20_phase_longread.sh`)
- `10`, `20`, `40`, `70`, `80`, `90` 사이의 빈 번호가 무엇을 뜻하는지 아무 데도 없었다

이제 **번호 = 파일 접두사 = 실행 순서**다.

## 매핑

| 이전 | 현재 | 단계 |
|---|---|---|
| `04_prepare_analysis_bam.sh` | `01_prepare_analysis_bam.sh` | WGS |
| `05_preflight.sh` | `02_preflight.sh` | WGS |
| `10_prepare_vcf.sh` | `03_prepare_vcf.sh` | WGS |
| `20_call_longread_clair3.sh` | `04_call_longread_clair3.sh` | WGS |
| `23_clair3_qc.sh` | `05_clair3_qc.sh` | WGS |
| `21_haplotag.sh` | `06_haplotag.sh` | WGS |
| `25_denovo_assembly.sh` | `07_denovo_assembly.sh` | WGS (de novo 분기) |
| `26_anchor_assembly.sh` | `08_anchor_assembly.sh` | WGS (de novo 분기) |
| `40_personal_fasta.sh` | `09_personal_fasta.sh` | WGS |
| `41_windows.sh` | `10_windows.sh` | WGS |
| `22_realign_wgs.sh` | `11_realign_wgs.sh` | WGS |
| `70_compact_zarr.sh` | `12_compact_zarr.sh` | WGS |
| `30_lift_annotation.py` | `13_lift_annotation.sh` (`.py` backend) | RNA 독립 단계 |
| `31_realign_rna_star.sh` | `14_realign_rna_star.sh` | RNA |
| `32_coverage_zarr.py` | `15_coverage_zarr.sh` (`.py` backend) | RNA 독립 단계 |
| `33_allelic_counts.py` | `16_allelic_counts.sh` (`.py` backend) | RNA 독립 단계 |
| `36_cnv_filter.py` | `17_cnv_filter.sh` (`.py` backend) | RNA 독립 단계 |
| `35_ase_qc.py` | `18_ase_qc.sh` (`.py` backend) | RNA 독립 단계 |
| `34_stamp_rna_metadata.py` | `19_stamp_rna_metadata.sh` (`.py` backend) | RNA 독립 단계 |
| `80_cleanup_intermediates.sh` | `20_cleanup_intermediates.sh` | 공용 |
| `90_make_catalog.py` | `21_make_catalog.py` | 유틸(수동 실행) |
| `20_phase_longread.sh` | `legacy_phase_longread.sh` | **번호 제거** — 파이프라인이 호출하지 않는 legacy다 |

번호 순서가 바뀐 곳은 셋이다: `22`(재정렬)가 `09`·`10` 뒤라는 사실이 이제 번호에 드러나고,
Clair3 QC(`23`→`05`)가 haplotag(`21`→`06`) 앞이라는 것도 드러난다.
RNA에서는 CNV 필터(`36`→`17`)가 ASE QC(`35`→`18`) 앞에 온다.

## ⚠ 깨지는 것

### `WGS_STOP_AFTER` 값이 전부 바뀐다

이전 값을 그대로 쓰면 **다른 단계에서 멈추거나 아예 멈추지 않는다.**

| 하고 싶은 것 | 이전 | 현재 |
|---|---|---|
| 입력 표준화까지 | `04` | `01` |
| preflight까지 | `05` | `02` |
| VCF 준비까지 | `10` | `03` |
| Clair3까지 | `20` | `04` |
| haplotag까지 | `21` | `06` |
| 개인 FASTA까지 | `40` | `09` |
| 재정렬까지 | `22` | `11` |
| compact까지 | `70` | `12` |
| cleanup까지 | `80` | `20` |

특히 위험한 값: **`20`은 이전에 Clair3였지만 지금은 cleanup(맨 끝)이다.** 옛 스크립트에
`WGS_STOP_AFTER=20`이 남아 있으면 멈추길 기대한 지점이 아니라 끝까지 돈다.

### work 디렉터리 로그 파일명

`logs/31.*` → `logs/14.*`, coverage는 `logs/15.cov.HAP*.log`, `logs/33.log` → `logs/16.log` 등이다. 옛 실행의 로그는 옛 이름 그대로 남는다.

2026-09-17 당시에는 RNA 13–19가 독립 wrapper였다. 이후 v2 RNA 단계가 13–26으로 확장·재배치되었으므로 현재 실행 순서는 문서 앞부분의 v2 표와 `run_rna_sample.sh`를 따른다.

### 바뀌지 않는 것

- 산출물 경로·파일명 (`individual_fasta/`, `compact/`, `allelic_counts.tsv` …)
- 설정 변수 이름 (`config_v2.sh`)
- `run_wgs_sample.sh` / `run_rna_sample.sh` / `run_pipeline.sh` 진입점
- `HANDOFF.md`의 2026-08 실행 로그 원문 — 기록이므로 옛 이름을 그대로 뒀다
