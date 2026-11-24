# process_code/v2 CPU·I/O 병렬 실행 정책

갱신일: 260928

## 목적

72 logical CPU 서버에서 검증 강도, read 선택, 정렬 파라미터와 최종 schema를 바꾸지 않고 전체 처리 시간을 줄인다. sample 수만 무작정 늘리지 않고 단계의 병렬 특성에 따라 동시 sample 수와 sample 내부 thread를 함께 제한한다.

## 기본 실행 정책

- 일반 `wgs-list`: 동시 WGS 2개, sample당 CPU 예산 36
- shard/대량 WGS: `wgs-phased-list`를 사용해 01–10은 4개, 11은 4개, 12–20은 2개 sample 동시 실행
- RNA: 동시 12개, sample당 CPU 예산 6
- 다운로드: shard 내부 2개 파일 병렬, shard queue에서는 현재 shard 계산과 다음 shard 다운로드를 중첩
- GPU: 기본 사용 안 함

```bash
# 대량 WGS 권장
./run_pipeline.sh wgs-phased-list /path/to/sample_list.txt all

# 단일 흐름 호환 실행
./run_pipeline.sh wgs-list /path/to/sample_list.txt all

# shard 목록의 다운로드/계산 중첩
./run_shard_queue.sh /path/to/shard_ids.txt \
  --data-root /data/mfm/input --output-root /data/mfm/output
```

## Stage 10: Illumina personal-HAP 재정렬

`all_to_both` routing은 indexed BAM/CRAM을 source header의 contig 순서로 나누어 읽고, 최대 `WGS_ROUTE_JOBS=16` process로 병렬 처리한다. 각 worker는 기존과 똑같이 secondary/supplementary를 제외하고, source duplicate를 `OD`, 원래 unmapped read를 `XU`로 운반한다. part는 원래 contig 순서로 결합한 뒤 기존 `samtools collate`와 pair QC를 그대로 거친다. HP 기반 `hp_routed_with_unknown`과 fragment 판정이 필요한 `hap_specific`은 의미가 바뀌지 않도록 기존 직렬 routing을 사용한다.

72코어 phased 기본의 Stage 10은 4 sample × 18 CPU 예산이다.

- routing: sample당 16 worker, 전체 최대 64 worker
- HAP alignment: sample당 HAP 2개 × (BWA-MEM2 7 + sort 2) = 18 thread
- 전체: 4 × 18 = 72 thread

개인 HAP FASTA는 sample마다 고유하므로 BWA-MEM2 index도 기본 `ephemeral`로 sample scratch에 생성하고 두 HAP 성공 후 제거한다. 수천 개의 개인 index가 전역 cache에 누적되는 것을 막는다. 반복 실험에서만 `STAGE11_BWA_INDEX_MODE=shared`를 명시한다.

`--hap-jobs`가 CLI config에 실제 연결되어 HAP1/HAP2 병렬 실행이 적용된다. 이전에는 인자가 전달돼도 config 생성에서 누락되어 기본값 1로 실행되던 결함이 있었다.
`config_v2.sh`는 첫 entrypoint의 auto/명시 요청을 `MFM_*_REQUESTED` provenance로 넘긴다. 따라서 run_pipeline→phased→runner→sample에서 config를 여러 번 읽어도 부모 fallback을 사용자 명시값으로 오인하지 않고, child의 18/36 CPU 예산으로 다시 계산한다.

재현을 위해 다음 변수를 명시할 수 있다.

```bash
export WGS_ROUTE_JOBS=16   # auto 또는 양의 정수
export WGS_HAP_JOBS=2      # auto|1|2
export WGS_AUTO_THREADS=1
```

## local scratch

공유 파일시스템의 임시 I/O를 줄이려면 노드 로컬 경로를 지정한다.

```bash
export MFM_SCRATCH_ROOT=/local_scratch/mfm
# 필요하면 단계별 override
export MFM_STAGE11_SCRATCH_ROOT=/local_scratch/mfm/stage11
export MFM_STAGE12_SCRATCH_ROOT=/local_scratch/mfm/stage12
```

Stage 10 routed BAM·sort 임시 파일과 Stage 12 compact 임시 파일만 scratch에 두며, 최종 산출물과 QC는 `OUTPUT_ROOT`에 기록한다. `MFM_STAGE11_SCRATCH_ROOT`라는 환경변수명은 구형 내부 Stage 11 호환을 위해 유지한다. sample output 절대경로의 SHA-256 prefix를 scratch key로 사용하므로 서로 다른 run의 같은 sample ID가 충돌하지 않는다.

## Stage 12: Compact Zarr와 full trace

`COMPACT_JOBS=auto`는 sample CPU 예산의 약수를 1–8 범위에서 고르고 `MFM_BLOSC_THREADS`를 나눈다. 36 CPU/sample이면 6 process × Blosc 6 thread이고 WGS 2개 합계는 72다.

full trace의 balanced 설정은 read metadata·packed sequence·CIGAR를 65,536-record 연속 범위로, high 설정은 1,048,576-record 범위로 prefetch한다. source SHA-256 multiset의 SQLite batch는 각각 8,192/262,144 rows다. 임시 DB만 `journal_mode=OFF`, `synchronous=OFF`, exclusive lock으로 만들고 완료 후 제거하며 최종 trace JSON과 진행상태 DB에는 적용하지 않는다.


## 자원 guard

runner는 CPU뿐 아니라 `MemAvailable`과 scratch 또는 output 파일시스템의 가용량으로 `LOCAL_JOBS`를 자동 하향한다.

- balanced/high WGS 기본 예약: 64/256 GiB RAM, 400 GiB 임시공간/sample
- balanced/high HAP 병렬 허용: 기본 64/256 GiB RAM, 300/400 GiB 임시공간/sample, CPU 16개/sample 이상 (운영체제·작업 특성에 따라 override 가능)
- balanced/high RNA 기본 예약: 48/64 GiB RAM, 80 GiB 임시공간/sample
- `/NHNHOME/korea` 가용량이 2 TiB 부근이면 `MFM_SCRATCH_ROOT`를 우선 사용

`LOCAL_JOBS=all|max`는 허용하지만 위 guard를 우회하지 않는다. core 수를 넘는 의도적 oversubscription은 기본 정책이 아니다.

## 실패 격리와 기록

한 sample의 preflight·정렬·검증 실패는 해당 sample의 `status.tsv`에 기록되고 다른 sample은 계속된다. phased 실행은 성공 sample만 다음 phase에 전달한다.

- sample 전체 자원: `run-status/.../logs/<sample>.resource.json`
- Stage 10 단계 시간: `qc/<sample>.stage11_timing.json`(구형 schema/파일명 호환)
- HAP별 자원: `logs/personal_realign/<sample>.HAP*.resource.json`
- Stage 12 contig 시간: `contig_timing.json`

## 모니터링

```bash
pgrep -af 'run_wgs_sample.sh|run_rna_sample.sh|bwa-mem2|minimap2|samtools sort'
uptime
free -h
df -h /NHNHOME/korea /tmp
```

메모리나 scratch 부족, 반복적인 tool failure, 높은 I/O wait와 완료율 저하가 나타나면 sample 동시 수를 낮춘다. 실행 시작 시 code snapshot을 고정하므로 원본 변경은 진행 중 batch에 소급 적용되지 않는다.


## 2026-09-29 재검증

- 현재 72-core 정책은 compact에서 `2 samples × 6 COMPACT_JOBS × 6 Blosc threads = 72`, RNA에서 `12 samples × 6 threads = 72`다.
- HAP 병렬 Stage 10은 sample당 `(14 align + 4 sort) × 2 HAP = 36`, 두 sample 합계 72코어를 넘지 않는다.
- `runner_local.sh`는 CPU 상한 뒤에 RAM과 scratch 상한을 함께 적용한다. 현재 서버처럼 `/NHNHOME/korea` 여유가 약 1.04TiB이면 WGS 400GiB/sample 기준 동시 sample은 최대 2로 자동 제한된다.
- 현 코드 chr22 full trace 실측은 HG00281 27.46초, HG00108 18.53초, NA20503 4.75초다. 세 sample의 QC와 manifest SHA-256은 2026-09-27 보존본과 각각 동일했다.

## 대용량 RAM 프로필 확정 (2026-09-29)

현재 서버 실측 가용 RAM은 약 1.9 TiB다. 기본 `auto`는 이를 `high`로 판정하고
약 20%(실측 391 GiB)를 OS/page cache/다른 작업에 남긴다. high의 WGS 예약은
256 GiB/sample이며, `/NHNHOME`에서는 기존 400 GiB scratch guard와 함께 적용된다.
scratch가 `/tmp` 같은 tmpfs이면 `sample RAM + sample scratch`를 합산한다.

3개 chr22 high-memory 재생성은 모든 Zarr array count·shape·dtype·logical bytes와
trace QC가 baseline과 동일했다. writer validation buffer를 4배 이상 키운 조합은
22.93→22.82초, 18.35→18.63초, 7.51→8.26초로 일관된 이득이 없어 채택하지 않았다.
반면 trace prefetch/SQLite 확대는 다음과 같아 high 기본으로 채택했다.

| sample | baseline trace | high trace | high max RSS |
|---|---:|---:|---:|
| HG00281 | 25.20초 | 19.54초 | 0.48 GiB |
| HG00108 | 15.90초 | 13.32초 | 0.36 GiB |
| NA20503 | 3.10초 | 3.04초 | 0.13 GiB |

HG00108 chr22 전체 14,217,688-read trace에서도 baseline과 high 결과 JSON이
byte-identical이었다. balanced 공유경로 임시 DB는 1,927.57초(32분 7.57초),
high batch/prefetch만 적용한 공유경로 실행은 1,161.54초(19분 21.54초),
high+임시 DB 최적화+`/tmp` scratch는 911.21초(15분 11.21초)였다. 최종값은
baseline 대비 **52.7% 단축(2.12배)** 됐고 max RSS는 각각 7.28/7.29/7.30 GiB였다.
baseline은 실행 중 물리 write가 163 GiB를
넘었지만 tmpfs 실행은 `/proc/<pid>/io` 기준 0이었다.

파일 open 한도는 650,000개로 충분하므로 open 수는 늘리지 않는다. 대형 RAM은
sort run 감소, trace 연속 읽기, SQLite transaction, OS page cache에 사용한다.

## Stage 12 trace process 병렬화 확정 (2026-09-29)

`TRACE_JOBS=auto`는 `COMPACT_JOBS`와 같은 worker 수를 사용한다. 각 작업은 `(contig,HAP)` 하나이며 단일 flat pool로 실행한다. 따라서 36 CPU/sample 기본은 `6 workers × 6 Blosc threads = 36`, WGS 두 sample은 최대 72 CPU다. 명시적 `TRACE_JOBS`도 `COMPACT_JOBS`로 cap한다.

3개 소형 chr22 pilot의 workers=2 단축률은 HG00281 21.9%, HG00108 46.9%, NA20503 49.6%였고 모든 JSON이 serial과 byte-identical이었다. HG00108 전체 14,217,688-read chr22는 911.21→459.48초로 49.6% 추가 단축됐다. 최초 balanced 1,927.57초 대비 총 76.2%(4.20배) 단축이다. worker별 peak RSS는 약 7.28 GiB, 관측 합산 RSS는 약 15 GiB였다.

worker 완료 순서와 무관하게 최종 report는 입력 contig와 HAP1/HAP2 순서로 쓴다. worker 하나의 예외는 구조화된 fail QC로 변환해 다른 HAP/contig 검증을 계속한다.

실제 `2 contigs × 2 HAP` smoke에서도 serial 29.36초가 workers=4에서 14.99초로 49.0% 단축됐고 contig QC 두 개와 manifest는 모두 byte-identical이었다.
