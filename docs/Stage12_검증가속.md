# WGS Stage12 독립 검증 가속

## 바꾼 부분과 바꾸지 않은 부분

Stage12는 compact 생성/내부 full 검증 후, source HAP CRAM과 저장된 Zarr를 독립적으로 비교한다.
HG00096 최초 실측은 생성·내부 검증 등38분15초 + 독립 trace1시간44분34초였다.
이번 개선은 trace의 계산 방식이며 alignment preset·read 필터·HAP mode·window·압축·배열 값은 바꾸지 않는다.
HG00096은 이미 cleanup됐고 source HAP CRAM이 없으므로 이 샘플의 전체 trace를 다시 실측했다고 하지 않는다.

| 대상 | 기존 | 변경 |
|---|---|---|
| read의 window 후보 | 매 read를 염색체 모든 window에 CIGAR 투영 | 정렬된 start/end 이진 탐색 후 겹칠 수 있는 window만 같은 projector에 전달 |
| depth | M/=/X 구간의 모든 염기에 매번 +1 | endpoints 최대8192개 batch, repeated index를 보존하는 add.at/subtract.at, 마지막 window별 누적합 |
| supplementary 검사 | segment SHA·Zarr 읽기·임시 DB를 또 반복 | 첫 pass에서 모든 segment를 검사하고 두 번째는 supplementary depth만 대조 |
| 임시 SQLite | rowid와 별도 PK B-tree | WITHOUT ROWID, worker별64MiB page cache |
| 진행 파악 | 작업이 끝난 뒤 전체 report | task 시작/종료, phase별 초, 병렬 task 완료 개수 즉시 flush 로그 |

trailing insertion은 reference end 위치를 다음 window 후보에 포함한다. CIGAR op/length 검사도 후보 선택 전에 수행한다.
손상되거나 단조 정렬이 아닌 window 목록은 기존 전수 탐색 경로로 가므로 QC를 숨기지 않는다.
depth는 D/N과 padding을 제외하고 primary/supplementary 정책을 그대로 사용한다.
uint32 modulo 누적 결과도 기존과 같다. read가 없는 window에는 불필요한 누적합을 하지 않는다.

검증기는 writer projection/depth 코드를 가져다 쓰지 않는다. 기존 schema·count·hash·audit·mismatch 예시 순서·PASS/FAIL 계약을 유지한다.
Zarr readset에 없는 base quality를 trace canonical hash가 검사한다고 주장하지 않는다.
SEQ/QUAL 보존은 Stage01 및 BAM→FASTQ/재정렬 입력 QC의 별도 검사다.

## 동일성·속도 확인

NA12878 HAP1 실제 chr22 재정렬 CRAM에서 coordinate-order record5000개를 임시 BAM으로 추출해 동일 입력을 비교했다.
51개1Mb window 중3854개 alignment가 유지됐고 supplementary51개였다.
source audit·canonical SHA·segment·depth mismatch 결과·SQLite count/예시가 백업 원본과 같다.

- source 투영: 원본2.418초 → 새0.538초, 이 작은 subset에서 약4.49배.
- primary depth/segment 대조:0.144초 →0.100초, 약1.44배.
- supplementary를 동일하게 segment까지 재검사하면0.063초 →0.066초로 동등 수준이다. 실제 새 경로는 중복 segment 검사를 하지 않아 별도로 측정한다.
- SQLite:10만 digest×source/output×read/segment=40만 insert + 의도적 mismatch1개.1.068초 →0.881초, 약1.21배. 임시 DB20,365,312bytes →10,264,576bytes.
- 최종 재측정도 모두 동일했다. source2.188→0.509초(4.30배), primary0.137→0.093초(1.47배), 중복 segment 검사를 제거한 supplementary0.059→0.024초(2.43배), SQLite1.027→0.866초(1.19배)였다. 이 수치는 해당 작은 component subset의 값이다.
- 최초 depth 구현은 빈 window에도 누적합을 해 sparse fixture에서 느렸고 peak RSS가 증가했다. 이를 제거한 최종 비교의 프로세스 전체 peak RSS는124,228KiB였다. 순차 실험의 누적 high-water 값이며 각각의 worker 독립 peak가 아니다.

추가 회귀는 window 경계 insertion/D/N/clip/P, dropped/overlap window, 잘못된 CIGAR, 반복 depth endpoints·uint32 wrap, duplicate/supplementary filter audit와 multiset multiplicity를 검사한다.
기존 fixture에는 read metadata·SEQ·CIGAR·segment·depth 변조 각각의 FAIL 및 원상 복원 PASS를 추가했다.
v2 반영 후 `tests/run_all.sh` 전체 **44 PASS / 0 FAIL**. 추가 가속 테스트6개와 기존 trace 테스트3개(각 corruption subtest 포함)도 모두 통과했다.

작은 warm-cache subset이며 전체 Stage12 또는 WGS가4.49배 빨라진다는 뜻은 아니다.
CRAM decode·Zarr unpack·SHA hashing·SQLite I/O 및 큰 염색체 straggler는 남아 있다.

## 자원 설정과 실행 중 적용

```bash
# 예: sample에14코어가 배정돼 있고, 메모리 실측이 충분할 때 다음 실행에서 사용.
MFM_SAMPLE_CPU_BUDGET=14 TRACE_JOBS=auto \
MFM_TRACE_CPU_RESERVE=1 MFM_TRACE_BLOSC_THREADS=1 \
./12_compact_zarr.sh SAMPLE
```

44개의 contig×HAP task 중 sample budget/사용자 max cap에 맞게 worker가 제한된다.
14코어·reserve1 조건에서는 최대13worker, 기본 reserve4에서는 최대10worker다.
기존 환경에서 `TRACE_JOBS=4`를 명시했다면 그대로4worker다. 실제 affinity와 부모의 총72코어 배분을 맞춰야 한다.
이 예시를 다른 sample와 CPU가 겹치게 독립 실행하라는 뜻은 아니다.

`MFM_TRACE_SQLITE_CACHE_MIB`는 기본64, worker별 cache다.13worker라면 최대 약832MiB cache에 source Python metadata/depth/decompression 메모리가 별도로 필요하다.
전체 DB를 `/dev/shm`에 무제한 넣거나 최종 출력 경로를 자동으로 바꾸지 않는다.
현재 scratch 기본값과 선택 옵션은 유지한다.

Stage12를 이미 실행 중인 프로세스는 로드한 이전 코드를 유지하고, 변경 반영 후 새로 뜨는 프로세스가 가속 코드를 사용한다.
기존 로그를 지우거나 실패를 PASS로 바꾸지 않는다. 다음 실측 로그로 source/SQLite/depth 병목과 CPU/RSS를 다시 판단한다.


## 2026-10-20 갱신

기본 multiset을 `auto`의 bounded 분할 정렬로 갱신함. SQLite는 임시 공간 부족 시 fallback 또는 명시 선택으로 보존함. read batch 4bit 해제·canonical digest 최적화, chunk 캐시64MiB, packed batch16MiB, scratch 사전 계산과 tmpfs RAM 계상을 추가함. SHA256 원값·중복 수·전수 depth 대조 기준을 유지함. 설정·재개·자원 한도는 [변경이력의 가속 정책](original_변경이력.md#rna-hap-병렬-처리독립-source-trace-가속-2026-10-20)에 통합함. 위의 10월1일 SQLite 측정은 과거 구성의 실측임.
