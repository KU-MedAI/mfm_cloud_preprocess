# 외부 서버 전달용 10TiB 계층화 전처리 안내

## 구성과 규모

이 전달본은 파이프라인, 공통 자료 목록, 3,212개 sample manifest와 shard 실행 목록을 하나의 폴더에서 관리한다.

- WGS–RNA pair: 449개
- WGS-only: 2,763개
- primary 입력: 80.890TiB
- URL 중복 제거 후 현재 크기를 아는 다운로드 자료: 약 80.935TiB
- 크기가 아직 동결되지 않은 CRAI/BAI와 1KGP panel 일부는 별도 표시되므로 실제 다운로드는 약 81TiB 이상이다.
- 현재 예상 최종 산출물: 307.309TiB
- 기본 분할: 37개
- 목표: shard당 8.5TiB
- hard cap: shard당 10TiB

최종 산출물 추정은 완료된 Illumina Compact Zarr, 기존 long-read 관측치, 입력 WGS/RNA 바이트를 사용한다. HiFi 전구간 및 RNA 전체 파일럿이 끝나기 전까지는 provisional 값이다.

## 폴더

```text
mfm_preprocess_share_260928_v4/
├── README.md
├── run_shard.sh
├── run_shard_queue.sh
├── pipeline/
├── catalog/
│   ├── shard_index.tsv
│   ├── sample_assignment.tsv
│   ├── stratification_qc.tsv
│   ├── runtime_estimate.tsv
│   └── shard_summary.json
├── shared/
│   ├── download_plan.tsv
│   └── manifests/
├── shards/
│   └── shard_001 ... shard_037/
├── data/
└── output/
```

공통 reference와 population VCF는 `shared/`에서 한 번만 관리한다. 각 shard의 `download_plan.tsv`에는 해당 shard sample 전용 WGS, RNA, index와 sample VCF만 기록한다.

## 72코어 시간 추정

다운로드 시간은 제외하고 72코어 서버 한 대를 계속 사용하는 조건이다. 기존 222.2일은 Stage 11 직렬 source scan/HAP 실행과 단계 공통 sample 병렬도를 사용한 보수적 기준선이다. 260928의 4→4→2 scheduler, contig routing, HAP CLI 연결, Stage 12 prefetch를 적용한 provisional 계획치는 다음 처리율로 계산한다.

| 구분 | 빠른 가정 | 중앙 가정 | 보수적 가정 |
|---|---:|---:|---:|
| Illumina WGS 처리율 | 30 samples/day | 24 samples/day | 20 samples/day |
| long-read WGS 처리율 | 18 samples/day | 12 samples/day | 10 samples/day |
| RNA 추가 처리율 | 192 samples/day | 144 samples/day | 100 samples/day |

| 구성 | 샘플 수 | 빠른 경우 | 중앙 추정 | 보수적 경우 |
|---|---:|---:|---:|---:|
| 1KGP Illumina WGS | 2,918 | 97.3일 | 121.6일 | 145.9일 |
| long-read WGS | 294 | 16.3일 | 24.5일 | 29.4일 |
| RNA 추가 단계 | 449 | 2.3일 | 3.1일 | 4.5일 |
| 전체(단일 서버 순차 합산) | 3,212 | 115.9일 | 149.2일 | 179.8일 |
| shard 평균(37개) | 약 87개 | 3.1일 | 4.0일 | 4.9일 |

따라서 현 시점 운영 계획은 222일이 아니라 **약 116–180일, 중앙 약 149일**이다. 서버 4대를 동일하게 쓰고 shared reference 다운로드 병목이 없으면 단순 계산상 약 29–45일이다. 다만 이 값은 다음 근거와 한계를 가진다.

- HG00373 과거 산출물 timestamp에서 personal realignment가 HAP1/HAP2 직렬로 약 4시간 10분, compact가 약 1시간 55분이었다.
- 최종 HG00281 pilot에서 contig route 16 worker가 약 1,600% CPU를 사용해 chr1–22와 unmapped의 23개 part 준비를 약 7분 20초에 완료했다. 변경 전 source scan은 sample당 약 1코어로 68분 뒤에도 진행 중이었다.
- 실제 standard collate는 scratch 사용과 시간이 커서 최종 구현은 8 thread와 `samtools cat | samtools collate -l 1` 스트리밍을 사용한다.
- BWA-MEM2 전구간, 압축 collate 여부, HiFi 전구간의 최종 성공 pilot이 아직 완료되지 않았으므로 149일은 확정 SLA가 아니다.
- 각 shard 완료 후 `resource.json`과 stage timing으로 플랫폼·프로젝트별 처리율을 다시 적합해 남은 ETA를 갱신한다.

## 실행

공유 루트에서 다음 명령 하나로 shard를 다운로드하고 pair 우선으로 처리한다.

```bash
./run_shard.sh shard_001 \
  --data-root /data/mfm/input \
  --output-root /data/mfm/output
```

실행 순서는 다음과 같다.

1. shard sample 전용 파일 다운로드 및 크기/checksum 검사
2. WGS–RNA pair의 WGS를 01–10(4 samples) → 11(4 samples) → 12–20(2 samples) 단계형으로 처리
3. 같은 pair의 RNA를 12 samples 병렬 처리
4. WGS-only도 같은 4→4→2 단계형으로 처리
5. 실패 sample을 상태표에 남기고 성공 sample만 다음 단계로 전달

명령만 확인하려면 다음을 사용한다.

```bash
./run_shard.sh shard_001 --mode preprocess --dry-run
```

다운로드만 하거나 검증만 할 수도 있다.

```bash
./run_shard.sh shard_001 --mode download
./run_shard.sh shard_001 --mode verify
```

여러 shard는 현재 shard 검증·전처리와 다음 shard의 2-way 다운로드를 중첩한다.

```bash
./run_shard_queue.sh /path/to/shard_ids.txt \
  --data-root /data/mfm/input --output-root /data/mfm/output
```

## 결과 동등성을 유지하는 속도 개선

적용된 개선은 정렬·필터·압축 파라미터를 바꾸지 않는다.

- Stage 12 trace 검증의 read metadata, sequence와 CIGAR를 65,536 record 연속 범위로 prefetch
- segment metadata를 window별 연속 범위로 읽어 per-scalar Zarr 접근 제거
- SHA-256 multiset SQLite insert는 기존 8,192개 `executemany` transaction 유지
- WGS는 01–10과 Stage 11에서 4 samples × sample당 18코어, Stage 12–20에서 2 samples × sample당 36코어로 자동 제한한다. RNA는 12 samples × sample당 6코어다
- pair WGS와 RNA를 분리 실행해 STAR 단계에서 72코어 활용률을 높임
- 개인 HAP BWA-MEM2 index는 sample scratch에서 임시 생성하고 성공 후 제거(shared는 반복 실험 opt-in)
- Illumina all_to_both source routing을 indexed contig별 최대 16 worker로 실행하고 source contig 순서로 재결합
- `--hap-jobs` CLI 연결 누락을 수정해 HAP1/HAP2 병렬 실행을 실제 적용
- Stage 11/12 임시 파일을 `MFM_SCRATCH_ROOT`의 local scratch에 둘 수 있음
- WGS 단계별 4→4→2 sample scheduler와 다음 shard 선다운로드 queue 적용

실제 chr22 재검증 결과는 다음과 같다.

| sample | 기존 trace | prefetch 적용 | 속도 향상 | 결과 |
|---|---:|---:|---:|---|
| HG00281 | 3,227.69초 | 26.01초 | 124.1배 | QC/manifest byte-identical |
| HG00108 | 2,098.81초 | 17.16초 | 122.3배 | QC byte-identical |
| NA20503 | 33.91초 | 4.38초 | 7.7배 | QC/manifest byte-identical |

따라서 기존 구현에서는 per-scalar Zarr 접근이 큰 병목이었다. 남은 source CIGAR depth 재계산은 검증 강도를 유지하기 위해 생략하지 않는다.

## 검증 계약

- 3,212개 sample이 정확히 한 shard에만 있어야 한다.
- WGS–RNA pair 449개와 WGS-only 2,763개 합계가 보존되어야 한다.
- 각 shard의 pair와 WGS-only sample 수 차이는 전체 비율 기준 최대 1개다.
- 모든 shard의 예상 출력은 10TiB 이하여야 한다.
- Stage 12 최적화 전후 canonical SHA-256 multiset, CIGAR, depth, metadata와 변조 탐지 결과가 같아야 한다.
- 실제 출력이 예측을 초과하면 완료 shard는 고정하고 미실행 shard만 다시 생성한다.
