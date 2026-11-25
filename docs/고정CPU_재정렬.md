# 고정 16CPU의 HAP 재정렬 활용 개선

## 적용 내용

- 대상: CPU HiFi의 자동 자원 설정이며, 샘플이 실제 16개 CPU의 소유권과 동일한 affinity를 가진 경우임.
- 기존: HAP1·HAP2 순차 실행, 각 `minimap2 -t 10`, `samtools sort -@ 0`, 입력 해독 추가 스레드 1개임.
- 변경: HAP1·HAP2 동시 실행, 각 `minimap2 -t 8`, `samtools sort -@ 1`, 입력 해독 추가 스레드 1개임. 두 HAP이 같은 샘플 소유 CPU 16개를 공유함.
- 도구와 부가 프로세스의 논리 CPU 계획 합계는 30임. 유휴·입출력 대기 구간을 겹치기 위한 스레드 구성이고, 예약 CPU와 실제 실행 가능한 CPU는 16개임. OS 스레드 개수나 실제 평균 사용 코어와 다른 값임.
- 입력 read·routing·preset·MAPQ/FLAG 기준·태그·HAP 좌표·CNV 계산은 그대로임. 스레드 구성과 동시 실행 순서만 변경함.

## 자원 보호와 호환

- `sample_cpu.owner`의 토큰, 살아 있는 소유 프로세스, 소유 CPU ID 개수와 현재 affinity의 일치를 검사함. 검증되지 않으면 기존 스레드 합계 기반 계획을 사용함.
- RAM은 동시에 활성인 두 HAP의 reference index·스레드·sort 버퍼를 모두 합산하고 기존 여유분을 더함. 부족하면 실행 전에 HAP 순차 처리로 복귀함.
- HAP 하나의 정렬이 끝나도 다른 HAP의 CPU/RAM은 반환하지 않음. 기존 HAP별 index·QC·완료 반환을 사용하고, 샘플의 CPU 소유권은 샘플 종료까지 유지함.
- 동적 CPU 모드, GPU, 다른 CPU 예산과 ONT/Illumina는 기존 계획을 유지함. 정렬·sort·HAP 수·입력 해독 스레드를 명시했거나 `MFM_WGS_CPU_BALANCE`를 명시하면 해당 설정을 우선함.
- 기본 고정 샘플 설정(`cpu_allocation_mode=fixed_sample`, `sample_cpus=16`)과 자동 스레드 설정이면 추가 인자 없이 적용함. 이미 생성한 정렬 script와 완료 결과는 자동으로 바꾸지 않음.
- RNA를 보류한 WGS 실행에서는 Stage01의 미사용 RNA 입력 검사를 제외함. RNA를 처리하는 실행과 `full` 모드에서는 기존 전체 해독·quickcheck를 유지함.

## 구현 위치

| 코드 | 역할 |
|---|---|
| `src/integrate_wgs_rna/realign/align_resources.py` | 소유권·affinity 확인 및 CPU/RAM 계획 결정함 |
| `src/integrate_wgs_rna/sample_cpu.py` | 샘플 CPU 소유권과 CPU ID를 유지하고 affinity를 상속함 |
| `src/integrate_wgs_rna/realign/align_memory.py` | 두 HAP의 최악 동시 RAM을 계산함 |
| `src/integrate_wgs_rna/realign/phases.py` | 실행 직전 계획을 장부에 예약함 |
| `src/integrate_wgs_rna/realign/hap_resources.py` | HAP별 후속 index/QC/반환을 관리함 |

`[REALIGN-PHASE]` 로그에서 `reserved_cpus=16`, `hap_jobs=2`, `threads=8`, `sort_threads=1`, `cpu_execution_policy=fixed_affinity_overlap`, `logical_pipeline_cpus=30`, `cpu_affinity_ids`를 확인함. 실제 CPU 사용량은 별도 계측값으로 판단함.

## 검증 범위와 시간 해석

- 2026-11-25의 동일 16CPU 부분 시험을 각 2회 수행한 결과임. 실행 순서를 뒤집어 반복했고, 평균 CPU는 CPU초 합계를 실행 경과시간 합계로 나눈 값임. RSS는 프로세스 트리의 관측 최댓값이며 GB는 10진 단위임.

| 구성 | 평균 경과시간 | 평균 실제 CPU | 관측 RSS 피크 |
|---|---:|---:|---:|
| HAP 순차, 각 `-t 10`, sort `-@ 0` | 90.00초 | 4.15코어 | 9.98GB |
| HAP 동시, 각 `-t 8`, sort `-@ 1` | 47.50초 | 8.10코어 | 19.03GB |
| HAP 동시, 각 `-t 10`, sort `-@ 1` | 47.85초 | 8.07코어 | 19.57GB |
| HAP 동시, 각 `-t 12`, sort `-@ 1` | 47.83초 | 8.08코어 | 19.43GB |

- 선택한 구성은 부분 입력에서 1.89배 빨랐음. HAP당 10/12스레드보다 메모리와 스레드가 적은 8스레드를 사용함. 이 표의 4.15/8.10코어는 부분 시험 값이며, 기존 HG00097 전량의 7.49코어와 구분함.
- 실제 HG00097 read의 부분 입력으로 BAM 해독부터 CRAM 출력까지 비교함. 두 작업 모두 같은 reference를 사용한 자원 배분 실험이며, 변이가 있는 실제 HAP 전량 재실측과 구분함.
- 충분한 batch가 생기도록 같은 실제 read를 반복 연결한 시험도 수행함. 반복 입력의 처리시간을 전체 donor의 완료시간으로 보고하지 않음.
- 같은 입력의 정렬 위치·FLAG·MAPQ·CIGAR·서열·품질·태그와 레코드 수를 정규화 SHA로 비교함. 공식 파이프라인의 소형 HAP 통합시험에서도 기존 출력과 QC 결과의 동등성을 확인함.
- HG00097의 기존 Stage10은 22,740초/CPU 170,315초로 평균 7.49CPU임. 같은 CPU 작업량을 유지한다고 가정하면 평균 12CPU일 때 3시간 57분, 14CPU일 때 3시간 23분, 16CPU일 때 2시간 57분임. 계산상 조건부 환산이며 개선 후 전량 실측값이 아님.
- 입력 읽기·단일 작업·압축 마무리는 항상 16CPU를 채울 수 없음. 전체 전처리의 다른 단계가 같은 비율로 빨라진다고 해석하지 않음.
- 동시 HAP은 순간 RAM과 I/O 요구량을 늘림. 샘플 4개 동시 실행의 처리량·RAM·I/O는 별도 전량 실측으로 판단해야 함. CPU affinity는 다른 프로젝트의 CPU 진입을 막는 독점 예약 기능을 제공하지 않음.

- RNA 보류 pair는 Stage01과02에서 미사용 RNA를 검사하지 않으며,02 계약에 `WGS_RNA_PAIR_RNA_DEFERRED`/`DEFERRED`를 기록함. 일반 RNA/full 입력 검증은 유지함. 실제22염색체 BAM/VCF 회귀로 검사함.
