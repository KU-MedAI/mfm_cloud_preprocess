# 파일럿 5개 오류 복구 가이드

2026-10-01 v2 수정 기록. 원본 데이터·기존 reference·완료된 WGS는 보존한다.

## Stage 06 reference와 donor 계약

WhatsHap은 BAM에 존재하는 전체 염색체와 미정렬 read를 기본 출력에 보존한다. Stage 02·06은 공통 `sample_alignment_reference()`로 현재 BAM을 실제 정렬한 reference를 선택한다. Stage 01이 raw uBAM/FASTQ를 새로 정렬했으면 analysis reference, 기존 BAM/CRAM 정렬을 보존해 변환했으면 source decode reference를 사용한다. `reused_canonical_bam`도 보존된 `alignment_policy.mode`로 구분한다. 개인 HAP 생성·WGS window는 기존 sequence reference를 계속 사용한다.

NA20503처럼 전체 contig가 남은 BAM에 autosome FASTA만 넘기면 chrX 처리 중 실패한다. 반대로 HG00096처럼 analysis reference에 새로 정렬한 BAM을 서로 다른 원본 decode FASTA와 비교하면 M5 false failure가 생긴다. 양쪽 오류를 함께 수정하고, 실제 reference의 22개 M5 검사는 그대로 유지했다. 서로 다른 source·analysis FASTA를 쓰는 raw BAM의 fixture도 추가했다.

`haplotag_contract.py`는 실행 전에 BAM 전체 contig의 reference 존재·길이를 확인하고, VCF가 단일 canonical sample인지 검사한다. 모든 BAM read-group의 SM은 하나의 donor여야 한다. Manifest에 `source_sample_id`로 명시된 다른 이름(HG001→NA12878)이 정확히 일치할 때만 `--ignore-read-groups`를 추가한다. BAM header는 원래 이름을 유지하고 판단 근거는 haplotag QC의 `input_contract`에 남긴다. 여러 donor 또는 알 수 없는 이름은 계속 실패한다.

`--regions`는 범위 밖·미정렬 read를 없앨 수 있어 이번 복구에 사용하지 않는다. 설치된 WhatsHap의 `normalize_user_regions`, `include_unmapped = regions is None` 및 출력 루프를 확인했고, 작은 실제 정렬 시험에서 chr1·chrX·unmapped read의 이름·서열·quality·flag가 모두 보존됨을 확인했다.

## Stage 06 입력폴더 기본값

Illumina short-read와 public VCF long-read는 준비된 표준 BAM을 `--alignment`로 명시한다. 입력폴더는 manifest의 `raw_root`, 환경변수 `RAW_ROOT`, 준비된 BAM의 부모폴더 순서로 선택한다. 예를 들어 `prepared_input/SAMPLE.analysis.bam`의 부모폴더는 `prepared_input`이다. 따라서 별도 샘플 실행기에서 `RAW_ROOT`를 지정하지 않아도 공용 Stage 06이 실행될 수 있다.

short-read는 BAM 파일 존재를 먼저 확인한 뒤 기본 입력폴더를 결정한다. 폴더 기본값을 정하는 것은 파일을 생성하거나 BAM을 다시 정렬하는 작업이 아니다. BAM·VCF가 없으면 중단하고, reference·donor·HP/PS 및 partial 출력 보존 계약은 그대로 적용한다. 실제 파일은 명시된 `--alignment`에서 읽는다.

## partial 출력 보존

`HAPLOTAG_REBUILD_PARTIAL=1`일 때만 partial 출력을 `haplotagged/partial_backup/<sample>.<timestamp>/`로 이동한다. BAM·index·summary·metadata뿐 아니라 `.part`·taglist·QC도 보존한다. BAM/index/summary/metadata/taglist/QC가 모두 존재할 때만 완성 출력으로 재사용한다. 원본 WGS는 이동·삭제하지 않는다.

## HG00281의 정확한 CRAM reference

기존 12GiB FASTA는 chr11 일부가 불일치하고 `verified/`의 3.1GiB FASTA는 chr2의 전체 M5가 불일치했다. 기존 chr2는 CRAM M5와 일치했다. 두 파일에서 CRAM header의 3,366개 @SQ M5·길이가 정확히 맞는 전체 contig를 선택하여 `_reference_cache/HG00281_cram_m5.fa`를 생성했다. chr2만 기존 FASTA에서 선택했고 나머지는 검증본에서 선택했다. 임의 base 변경이나 MD5 검사 해제는 하지 않았다.

모든 @SQ를 검증했고 chr1–22의 전체 M5를 다시 대조했다. source/analysis reference는 파일럿 Manifest의 HG00281 행에만 이 자산을 연결했다. 과거 두 FASTA와 Manifest 백업은 보존했다. 검증 보고서는 출력 루트 `_run_status/261001_pilot5_recovery/recovery_contracts.json`에 있다.

## RNA-only 재개

HG00096 WGS compact와 chr1–22 trace는 완료되어 다시 생성하지 않는다. 빠져 있던 RNA additional FASTA는 기존 `15_star_index.py`의 `ensure_addnl`을 재사용해 173개 contig로 준비했다. source FAI와 additional FAI도 존재한다. 재개 시 Stage 02 계약을 다시 확인하고 RNA Stage 13부터 시작한다. RNA 완료와 sample package 검증 전에 전체 완료로 집계하지 않는다.

RNA-only 재개에는 기존 WGS 공통 preflight의 500GiB workspace 값 대신 명시적 128GiB RNA 예산을 전달한다. 전체 5개 WGS+RNA 재개는 2,048GiB 별도 guard를 유지한다. reference·sample·pair·quality 등의 과학적 검증 기준은 변경하지 않았다.

## 5개 동시 실행 정책

- 14 CPU × 5 worker = 70 CPU, 전체 예산 72 CPU. OS CPU affinity도 worker마다 겹치지 않게 설정한다.
- WGS HAP 내부는 직렬로 하여 sample 간 병렬과 중첩 과점유하지 않게 한다. aligner 10 + sort 2, compact 2 worker, trace 4 worker 등 sample 예산 안에서 제한한다.
- 상태·로그·결과 JSON을 sample마다 분리한다. 한 sample의 실패가 나머지 worker를 종료하지 않는다.
- 기존 work RUN_ID는 재사용하고 새 상태 기록은 별도 run label에 둔다.
- 모든 5개를 시작하기 전 여유 공간 2,048GiB 운영 예산을 검사한다. 이는 남은 4개 final WGS·haplotag BAM·sort·RNA 공간을 위한 안전 예산이지 실제 최종 크기 측정값은 아니다. 현재 약 408GiB여서 5개 동시 production은 차단됐다.
- HG00096 RNA-only는 128GiB 여유 공간 guard를 별도로 적용한다. 기존 공유 출력 경로는 바꾸지 않는다.

5개 실행용 스크립트는 `codex_archive/plan/261001/pilot5_동시재개.py`와 `.sh`다. 저장공간이 확보된 후 새 `--run-label`을 지정해야 하며 기존 실행 기록을 덮어쓰지 않는다.

## 테스트

기존 전체 회귀 시험 40개 스크립트가 모두 통과했다. 새 donor/reference/partial 시험은 10개이고 실제 WhatsHap read 보존 시험도 포함한다. donor 불일치·여러 sample·누락 reference·길이 불일치·partial 미보존을 계속 차단하는지 확인했다. 이 결과는 5개 전체 production 완료를 뜻하지 않는다.

## 14:02 업데이트 — window 계약 복구와 RNA 18 재개

HG00096의 구형 `windows.tsv`에는 `raw_start/actual_length`만 있어 RNA 18 reader가 실패했다. 공통 `src/integrate_wgs_rna/window_contract.py`에서 `personal_start=raw_start`, `personal_end=raw_start+actual_length`를 검증·보충한다. 두 표현이 함께 있으면서 다르면 중단한다. Stage 11의 캐시 재사용도 같은 검사 후에만 허용한다. 개인 좌표를 GRCh38 좌표로 복사하지 않으며 없는 REF 투영은 -1로 남긴다.

원본은 archive 백업과 sample window 폴더의 고유 `.legacy.*` 파일에 보존했다. 기존 5,361개 행·모든 기존 열 값·kept 선택은 그대로다. 22개 염색체 × 양 HAP, 총 44개 kept 구간 목록이 compact와 일치했다. compact index는 retained 행을 재번호화하므로 전역 window index와 번호가 다를 수 있다. Stage 23에서는 **좌표·길이·순서**를 대조하고 RNA window에서는 전역 key와 실제 좌표 배열·누적 offset까지 검사한다. WGS Zarr와 source trace는 수정하지 않았다.

전체 회귀 검사 42개 script PASS / 0 FAIL. 구형/신형 입력의 depth·개인 좌표·key 동일성, REF 미정 sentinel, 좌표 충돌·중복·overlap·변조 차단과 종료 코드/실패 로그 보존을 검사했다. 보고서는 sample의 `qc/window_contract_upgrade.json`과 archive의 `window로그_회귀검사.log`다.

HG00096은 `--rna-start-at 18 --run-label 261001_HG00096_rna_resume_r3`로 재개하여 Stage 02를 통과했고 13–17을 건너뛰었다. Stage 18 전장 coverage부터 재계산하므로 기존 HAP1 coverage도 다시 계산한다. 전체 완료로 간주하지 않으며 Stage 23 PASS 후에만 기존 cleanup 정책을 적용한다. 나머지 4개는 여전히 용량 부족으로 production 대기다.

## 실행 즉시 로그와 실패 이력 보존

예전부터 도구 로그·runtime·provenance는 기록했지만 work 로그의 final 수집은 Stage 23에 도달해야 수행됐다. 이제 공통 runner가 실행 시작부터 stdout/stderr를 최종 sample에 직접 기록한다. 원본 도구 로그도 기존대로 유지하고 Stage 23에서 추가 수집한다.

```text
<SAMPLE>/logs/stage_logs/<WGS|RNA>/<stage>.<고유번호>.log
<SAMPLE>/metadata/execution_provenance.jsonl
<SAMPLE>/metadata/stage_runtime.tsv
<SAMPLE>/metadata/stage_failures.jsonl
<SAMPLE>/metadata/recovery_history.jsonl
<SAMPLE>/logs/recovery/<원본경로별 고유번호>.log
<SAMPLE>/logs/runner/<run-label>/run.log
```

- stage 로그: 실행 시도마다 새 파일. stdout·stderr 동시 기록, 기존 파일 덮어쓰기 없음.
- provenance: 실제 argv·도구 버전·시각·상태·종료 코드와 해당 로그 경로 연결. runtime TSV는 실패 행의 reason에도 로그 경로를 남긴다.
- 실패 이력: bounded log tail·오류 문구·종료 코드·명령·정확한 로그 경로를 JSONL에 누적한다. 자동 추출은 `diagnosis_status=unclassified`이며 확인한 원인과 구별한다. 로그 기록 실패도 정상 성공으로 처리하지 않는다.
- 복구 이력: 기존 파일럿 5개 원인·수정·production 대기/재시도 결과를 기록한다. 테스트 PASS를 production 완료로 간주하지 않는다. 과거 로그 사본은 원본 위치와 SHA-256를 연결해 보존한다.
- 파일럿 controller는 stage runner 진입 이전 preflight 실패도 외부 `run.log` 사본으로 보존하고 시작/종료를 recovery history에 추가한다. 현재 실행 중인 외부 로그는 `_run_status/<run-label>/<sample>/run.log`이며 final runner 사본은 종료 시 보존한다.

## 14:32 이후 사용자 승인 — 고정 용량 시작 기준 폐지

위의 2,048GiB/500GiB/128GiB 설명은 과거 실행 정책이다. 사용자가 고정 예산 없이 실제 진행하도록 승인하여 현재 기본은 **가용 공간 관찰만 수행**한다. controller는 `minimum_free_gib=null`, `disk_policy=observe_only`로 실제 여유를 기록하며 작은 여유 공간 자체를 이유로 시작을 차단하지 않는다. Stage 02도 기본 `PREFLIGHT_MIN_FREE_GIB=0`, `policy=observe_only`다. 사용자가 명시적으로 최소값을 전달한 경우에만 해당 값을 적용한다.

파일시스템 조회 오류·실제 ENOSPC/I/O 실패·reference/SEQ/QUAL/donor/pair/최종 package QC는 그대로 유지한다. 디스크가 실제로 부족할 수 있으므로 이를 production 성공 보장으로 해석하지 않는다. 입력·성공 결과를 임의로 삭제하거나 출력 경로를 바꾸지 않는다. 나머지 4개는 CPU 14–69를 쓰며 기존 HG00096 CPU 0–13과 겹치지 않는다.

## 15시 이후 승인 — RNA 실행 효율 개선

- Stage 21 `21_cnv_filter.py`는 HAP×contig 44개 작업을 `rna_cnv_depth.py`에 맡긴다. bedcov의 필터·D/N 처리·구간은 그대로이고, 원래 BED 행 순서로 합산한 뒤 **HAP 전체** 유전자 중앙값을 사용한다. 다른 gene의 겹치는 exon을 합치거나 RNA에서 발현된 gene만 남기지 않는다.
- Stage 18 전장 `18_coverage_zarr.py`는 `rna_coverage.py`의 독립 contig process로 실행한다. Parent가 모든 dataset을 먼저 만들고 worker는 자신에게 배정된 contig 배열만 쓴다. 원래 raw uint32 depth·shape·chunk·zstd5·정규화 상수를 유지한다. 전체 성공 후에만 완료 attrs를 기록한다.
- Stage 18 window는 `--coverage-source-zarr`로 같은 전장 coverage에서 window를 잘라 저장한다. BAM/FAI path·size·mtime·inode, sample/HAP, 필터 인자·좌표계·shape를 대조하고 불일치는 출력 생성 전에 중단한다. genome→window 축약 시 모델의 좌표·행 순서·offset·raw 값은 바뀌지 않는다. 옵션을 생략하면 기존 독립 BED 계산 경로가 남아 있다.
- 자동 CPU 정책: 실제 affinity와 `MFM_SAMPLE_CPU_BUDGET`/`MFM_CPU_BUDGET`의 최소값을 사용한다. 조정용 1코어를 남기고 coverage는 Python+samtools 2코어/job, CNV는 1코어/job으로 계산한다. 샘플 14코어에서는 CNV 최대13, coverage 최대6 worker이며 Blosc 내부 worker는1이다. `RNA_CNV_JOBS=1`, `RNA_COVERAGE_JOBS=1`로 직렬 비교가 가능하고 지나치게 큰 명시값도 예산으로 제한된다.
- task 시작/완료·구간 수·경과 초·완료/전체 수가 sample stage 로그에 즉시 기록된다. 실패한 CNV task가 나오면 이 호출의 나머지 bedcov도 종료한다. 성공 TSV 없이 RNA22로 진행하지 않는다.
- HG00096은 성공한 RNA13–20을 다시 만들지 않고 RNA21부터 새 run label로 재개한다. 원래 직렬 실행 로그/실패/복구 이력은 보존한다. 현재 코드 변경은 v2뿐이며 v3는 원본 그대로다.

### 동일성 검사와 측정 범위

새 7개 회귀 케이스는 BAM/CRAM, exon 중복·겹침·접합·cross-contig gene, 짧은 gene/NA·0중앙값 실패, duplicate/secondary/supplementary/QC/MAPQ·D/N·mate overlap·무coverage contig, 독립 BED→window 대조, 입력 변조·외부 도구 실패·CPU 제한을 확인한다. 백업한 원본 Stage18의 raw 배열/기존 attrs와 Stage21 TSV도 동일하다.

실제 HG00096의 HAP1/HAP2 chr21/chr22 각 원본 exon 200행씩으로 기준선과 비교했다. 병합 후 합계236 BED행, 출력21 gene·999bytes의 SHA-256은 모두 `261d0e1044dc3e04c54f1863b9a573e4126709ce4477004ba62399b5bdc19b96`이다. 최초 원본8.759초, worker1 8.634초, worker4 2.914초였다. 작은 warm-cache subset 측정이지 42만 전체 구간·전체 RNA의 3배 개선 보장은 아니다.

HG00096 WGS01–12 성공 단계 합계는19,972초(**5시간32분52초**), 최초 시작부터 완료까지 wall-clock은8시간53분8초다. 후자는 오류 수정·재개 대기를 포함하고 별도 source-trace 재검사 및 RNA 후 package 검증은 이 합계와 구별한다.

### 속도 개선 테스트가 발견한 별도 정책 문제

2026-10-01 감사 당시 `rna_profile.py`/coverage attrs에는 duplicate 포함으로 적혀 있었지만 실제 `samtools depth -a -J -Q 255 -G 2820`는 **기본 DUP 제외를 유지**했다. `-G`는 제외 flag를 추가하며 기본 DUP 제외 해제에는 `-g 1024`가 필요하다. 설치된 samtools1.24 도움말·duplicate fixture·[공식 depth 문서](https://www.htslib.org/doc/samtools-depth.html)를 대조했다.

후속 `pipeline1_update` 통합에서 방향을 GTEx target과 동일하게 확정해 `-g 1024`를 추가했다. 따라서 Stage 18 coverage는 duplicate를 포함하고, 구형 duplicate-excluded coverage는 재사용하지 않는다. STAR index·정렬은 alignment-only fingerprint로 재사용하며 Stage 18–19만 다시 계산한다. ASE duplicate 제외는 유지한다.

## 추가 승인 — WGS Stage12 독립 검증 가속

HG00096은2026-10-01 15:57:28에98개 최종 package 검사 PASS 뒤 cleanup까지 완료했다.
추가 개선은 완료 출력을 다시 만드는 것이 아니라 v2의 이후 Stage12 실행에 적용한다.
window 후보 탐색·depth endpoints·SQLite 구조/cache·중복 supplementary segment 검증을 개선하고
read/segment hash, 필터, window 좌표, depth와 QC 실패 판정은 보존한다.
task/phase elapsed와 완료 개수는 기존 통합 로그에 남긴다.
자세한 동일성·성능 범위와 사용자 설정은 [Stage12_검증가속.md](Stage12_검증가속.md)에 기록한다.
