# original 변경 이력 참고본

## 2026-11-25 — 기본 pileup36채널

- 최신 기본: 공통12+HP1별9+HP2별9+공통BAF3+non-ref2+공통reverse1임. 신규 Short/Long 기본 출력은 pileup_feature임. readset full/sequence_only는 명시 선택함.
- read 필터·HAP/window·RNA/ASE/CNV 정책은 유지함. 단일 native 순회로 HP를 함께 집계하고 독립 HTSlib로 전수 대조함.
- 새 layout/schema v2를 사용하며 구버전12채널은 읽기 호환을 유지함. 신규 완료 판정은36 feature_set까지 확인함. [현행 상세](pileup_36채널.md)를 참조함. 아래 날짜별 절은 당시 이력임.


이 문서는 복사 시점 original README의 근거/과거 이력임. 클라우드의 실제 경로·기본 설정·실행법은 최상위 README와 docs/클라우드_실행.md를 우선함.

# MFM WGS·RNA 전처리 v2

## 현재 사용 안내 · 2026-11-23

- **처음 실행할 때:** [GUIDE.md](../GUIDE.md)의 모드 선택·JSON·dry-run/run·자원·정리 순서를 따름.
- **단계와 도구:** [PIPELINE.md](../PIPELINE.md), [단계별 입출력 상세](../docs/단계별_입출력_상세.md)를 확인함.
- **결과와 상태 경로:** [실행·출력 구조](../docs/실행과_출력구조.md), 외부 서버는 [PORTABLE_GUIDE.md](../PORTABLE_GUIDE.md)를 확인함.
- 코드의 WGS 기본 출력은 `full`임. `sequence_only`·`pileup_feature`는 JSON에서 명시 선택하며, WGS 모드와 RNA `full|coverage_only`는 별개임. `run_pipeline.sh run`의 내장 큐를 사용하며 별도 운영 제어기를 작성할 필요 없음.
- **WGS-only pileup은 개인 HAP 재정렬을 하지 않음.** reference 정렬·haplotag BAM에서 tensor를 계산하고 chain으로 HAP 서열 좌표에 투영함. `MFM_PILEUP_BENCHMARK_REALIGN=1`은 Stage10을 더하는 실측 전용 옵션임. 일반 처리에는0을 사용함.
- RNA `full`의 ASE/CNV가 필요한 pair는 WGS 근거 처리와 개인 HAP CRAM이 별도로 필요함. `coverage_only`는 HAP별1bp signal·QC만 만들며 ASE/CNV를 비적용으로 기록함.
- 현재 제공 자원 예시는 `fixed_sample`/샘플16CPU/최대4명/공유68CPU/RAM1.65TB·예약1.485TB임. 실제 더 작은 할당을 우선함. 동시수1만 지정해도 샘플 예산16이68로 늘어나지 않음. 배정량·실제 사용량·도구 스레드를 구분함.
- 모델 전용 cleanup은 최신 QC/trace 통과 후 새 중간 BAM/CRAM/FASTA/STAR index 등을 정리할 수 있음. HAP 서열 Zarr·window 주석·압축 chain·QC·metadata·로그를 보존함. 기존 원본을 보존할 실행은 [정리 설정 조합](../GUIDE.md#7-검증재개정리)을 확인함.
- 이 문서의 **이전 날짜별 절은 당시 구현·실험 이력**임. 과거1TB/900GB, 동적 대여, GPU 비교 예외, 당시 완료 인원·경로를 현재 기본값으로 사용하지 않음. 현재 사용법은 위 안내와 해당 실행의 `effective_config.json`·worker 장부를 기준으로 함.

<a id="reference-pileup-선택-모드-2026-11-17"></a>

## Reference pileup 선택 모드 (2026-11-23)

- 선택값: `pipeline_env.MFM_SHORT_READ_OUTPUT_MODE` 또는 `MFM_LONG_READ_OUTPUT_MODE`를 `pileup_feature`로 지정함. 기존 `full`·`sequence_only`와 기본값은 유지함. `acceleration.mode`·검증의 `--mode wgs|rna`와 다른 출력 형식 설정임.
- 보존: 기존 consensus HAP1/2 서열, chr1–22, 1,048,576bp window, seed42, 말단 N trim, 내부 N·blacklist 정책, padding·순서·좌표 주석을 유지함. 서열에 read pileup을 추가하며 최종 WGS readset은 저장하지 않음.
- 계산 기준: 기존 reference 준비 BAM 및 read 필터를 사용함. HAP에 WGS read를 다시 정렬한 값이 아니며 reference 관측을 기존 consensus chain으로 HAP window에 대응함. HAP1/2의 동일 reference 위치는 같은 공통 관측을 사용함. HP별·allele별 지지 채널로 해석하면 안 됨.
- 처리: 기존 WGS 준비·사전 검증·VCF·haplotag·HAP FASTA/chain은 유지함. WGS-only 및 RNA `coverage_only`의 기본 pileup에서는 WGS Stage10 HAP 재정렬이 비적용임. 실측 전용 `MFM_PILEUP_BENCHMARK_REALIGN=1`은 예외이며 일반 처리에는0을 사용함. RNA full의 ASE/CNV가 HAP CRAM을 요구하면 Stage10을 유지하며 해당 계산 정책을 변경하지 않음. `vcf_consensus`의 검증된 plus/plus chain만 지원함.
- RNA: STAR·HAP별 1bp coverage·annotation lift는 기존 경로를 사용함. `MFM_RNA_OUTPUT_MODE=coverage_only`는 기존 정책대로 ASE/CNV를 생략하는 별도 선택임. 이 WGS 출력 모드 자체가 RNA 정책을 바꾸지는 않음.
- read 필터: 플랫폼별 기존 필터와 `MIN_MAPQ`를 사용함. BQ cutoff·BAQ·depth cap을 추가하지 않음. 겹치는 paired read는 각 alignment로 셈. Long supplementary는 primary feature에 합산하지 않고 별도 depth로 기록함. missing QUAL/BQ255는 BQ 평균의 분모에서 제외함.

저장 위치는 `compact/hap_sequence_reference_pileup_v1/<chr>/`임. 기존 `hap_sequence.zarr` 및 window 원장과 함께 `pileup.zarr/HAP1|HAP2`에 아래 배열을 저장함.

| 배열 | 형식·의미 |
|---|---|
| `features` | float32, `[전체 padded position, 12]`임. 채널 순서는 아래 표와 manifest에 기록함 |
| `depth` | uint32, reference M/=/X 염기 관측수임. deletion span은 포함하지 않음 |
| `baseq_depth` | uint32, 실제 BQ가 존재하는 관측수임 |
| `supplementary_depth` | uint32, 허용된 supplementary의 별도 염기 관측수임 |
| `max_deletion_length` | uint32, anchor 뒤 read CIGAR deletion의 최대 길이 원값임. feature의 float32 표현과 별도로 정확한 정수를 보존함 |
| `mapping_kind` | uint8, 0=padding, 1=reference 대응, 2=HAP 삽입임. depth0과 비대응을 구분함 |
| `offsets`, `lengths` | 기존 window와 같은 위치·순서의 구간 인덱스임 |

| 순서 | 채널 | 정의 |
|---:|---|---|
| 0 | `log1p_depth` | `log(1 + depth)`임 |
| 1 | `mean_mapq` | 염기 관측 read의 MAPQ 평균임 |
| 2 | `mean_baseq` | 품질값이 있는 관측의 BQ 평균임. 관측 부재는 `baseq_depth=0`으로 구분함 |
| 3 | `strand_forward_ratio` | 염기 관측 중 forward 비율임 |
| 4–8 | `read_A/C/G/T/N_ratio` | A/C/G/T/N 관측수 각각을 base depth로 나눈 값임 |
| 9 | `insertion_ratio` | 왼쪽 anchor 뒤 insertion 사건수 / (anchor base depth + deletion span depth)임 |
| 10 | `deletion_after_event_ratio` | anchor 뒤 deletion 사건수 / 같은 분모임 |
| 11 | `deletion_after_max_len` | anchor 뒤 read CIGAR deletion의 최대 길이임. HAP↔reference chain의 deletion 길이가 아님 |

- HAP 삽입: 삽입 염기는 HAP 서열에 보존함. `mapping_kind=2`이며 왼쪽 reference anchor의 generic `insertion_ratio`만 전달함. 해당 HAP allele과의 일치율이나 삽입 염기별 BQ·base ratio·depth를 생성하지 않음. reference 왼쪽 anchor가 없으면 비율도 미관측0임. padding의 모든 값은0임.
- 정밀도: 비율·평균은 float32임. count와 deletion 길이 원값은 uint32임. 모든 서열·feature는 zstd level5로 압축함.
- 검증: reference/decode/HAP 생성 reference의 contig 서열 MD5·BAM SQ M5(있을 때)·chain 길이를 확인함. SQ M5가 없는 BAM은 기존 Stage02의 검증된 reference 근거에 의존하며 BAM만으로 정렬 reference를 증명할 수는 없음. source/policy 지문과 전 배열 SHA를 검사함. 생산 C++ CIGAR 순회와 독립 HTSlib pileup 집계로 저장 token 전체를 대조함. backend `auto`의 빌드/로드 불가 때만 기존 Python/aligned-pairs 경로로 대체하며 입력/계산 오류는 숨기지 않음. 최종 QC 통과 후에만 기존 cleanup을 수행함.
- 자원: sample CPU grant 안의 contig 병렬화와 bounded reference block cache를 사용함. worker 기본 RAM 예약은4GiB(`MFM_PILEUP_WORKER_RAM_BYTES`), 메모리 보호가 활성일 때 가상 주소공간 상한은8GiB(`MFM_PILEUP_WORKER_AS_BYTES`)임. 주소공간은 RAM 예약 이상이어야 하며 실제 RSS와 구분함. 극단적인 read 길이·밀도에서 실패하면 필터나 depth를 자르지 말고 RAM 예산을 재검토해야 함.
- 읽기: `WindowReader(sample_root).load(contig, hap, index)`에서 기존 DNA·RNA 결과와 `pileup_features`, `pileup_channels`, reference 대응/삽입/BQ 관측 mask, `max_deletion_length`, supplementary depth를 반환함. `wgs_depth_coordinate_system=reference_projected_to_hap`이며 full 모드의 개인 HAP 재정렬 depth와 구분함.
- 재개: 같은 입력·정책·완료 배열만 재사용함. 지문이 다르면 덮어쓰지 않고 실패함. 기존 완료본과 다른 `output_root`·`state_root`·`run_id`를 지정해야 함. 원본이 정리된 완료본에 feature를 자동 추가하지 않음.

기존 JSON을 복사하여 새 출력/상태 경로를 지정하고 다음 값을 `pipeline_env`에 병합함. pair의 기존 RNA 정책을 유지하려면 `MFM_RNA_OUTPUT_MODE`는 현재 설정을 그대로 사용함.

```json
"pipeline_env": {
  "MFM_SHORT_READ_OUTPUT_MODE": "pileup_feature",
  "MFM_LONG_READ_OUTPUT_MODE": "pileup_feature",
  "MFM_PILEUP_BENCHMARK_REALIGN": "0",
  "MFM_RNA_OUTPUT_MODE": "full",
  "MFM_SEQUENCE_DEFER_RNA": "0"
}
```

```bash
conda activate jj_mfm
bash run_pipeline.sh dry-run --selection /path/selection.tsv --config /path/pileup.json --executor local
bash run_pipeline.sh run --selection /path/selection.tsv --config /path/pileup.json --executor local
```

- 위 JSON은 RNA 기본 `full`을 유지함. ASE/CNV 없이 signal만 처리할 때에만 `MFM_RNA_OUTPUT_MODE=coverage_only`로 바꿈. RNA 보류와 동시에 사용하지 않음.

- 소스 위치: `src/integrate_wgs_rna/pileup_features/`의 `compute.py`는 계산, `projection.py`는 chain 대응, `export.py`는 저장·source 검증을 담당함. Stage11/12와 기존 검증·reader·cleanup에 연결되어 별도 연속 제어기는 필요하지 않음.
- 검증 범위: 합성 BAM의 indel·N·품질·padding·재개·중단·정리 및 RNA 결합을 검사함. 실제 HG00706 chr22의 HAP당 기존1Mb 창1개씩에서 기존 완료본의 창 좌표·DNA 코드 전체 일치와 독립 source 대조 PASS를 확인함. 두 창의 생성은32.8초, 독립 검증 포함268초·단일 프로세스 관측 최대 RSS617MB·compact17.25MB였음. 국소 최대 depth8,502인 부분 시험이며 전체 샘플 시간·용량의 실측값으로 해석하지 않음.

### 최근 전수 실측과 메모리 보호 보완

- 2026-11-20~21의 기본 pileup 전과정 실측은 동일 HiFi WGS-only4명, 고정16CPU×4/공유68, 다운로드 제외·중간자료 보존 조건임. 샘플별 약1시간56분~2시간11분, 최종 모델입력49.78~54.19GB를 기록함. HG01081은 실패/복구 실행시간을 합산하고 수동중단은 제외함. 모든 샘플/플랫폼의 처리시간 보장값은 아님.
- 과거 readset의 동적 공유CPU 기록과는 배분/필수 단계/정리 조건이 달라 저장 형식만의 가속률로 해석하지 않음. 강제 HAP 재정렬을 포함한 별도 실험 시간도 기본 pileup 시간과 구분함.
- HG01081 native 검증에서 가상 주소공간4GiB 제한의 실패와8GiB 통과를 확인해 worker 기본 보호를 RAM예약4GiB/주소공간8GiB로 보완함. 전역 RAM 부족과 구분하며 feature·필터·전수검증 범위를 바꾸지 않음.
- 아래 2026-11-18/19 기록의 worker2GiB·주소4GiB는 당시 설정임. 현재 값은 위 설명과 `pileup_features/export.py::_worker_memory_limits`를 따름.

## Short-read WGS 서열 전용 저장 (2026-11-02)

- 선택: JSON `pipeline_env.MFM_SHORT_READ_OUTPUT_MODE=sequence_only`로 실행함. 기본값은 `full`임. Short는 제공된 phased Illumina panel, Long은 `public`·`provided` VCF로 서열 전용 처리가 가능함.
- **WGS-only**: WGS BAM/CRAM을 다운로드하거나 읽지 않음. 동봉 공통 자산만 확인할 때는 원본 전송 동시수/디스크 대기와 분리하며 RAM 예약과 worker의 SHA 검증은 유지함. 개인 VCF 준비 → GRCh38 기반 HAP1·HAP2 FASTA/chain → 기존 window 생성 → 서열 Zarr → full/source QC → 중간 자료 정리 순서임.
- **Window 정책 보존**: chr1–22, 1,048,576bp, seed42의 결정적 시작점, 말단 연속 N trim을 유지함. 내부 N·blacklist 겹침을 보존하며 실제 구간이 전부 N인 창만 제외함. 기존 window 원장을 그대로 사용하며 N 비율·blacklist 비율/겹침·좌표 매핑을 보존함. 이 모드 때문에 window를 추가로 제외하지 않음.
- 저장: `compact/sequence_4bit_acgtnp__windowed/<chr>/hap_sequence.zarr`임. 기존 서열과 같은 4bit ACGTNP·zstd5 및 window padding/offset 형식을 사용함. `readset.zarr`와 `depth.zarr`는 생성하지 않으며 manifest/sample metadata에 `wgs_payload=sequence_only`로 명시함.
- **RNA pair**: RNA read·STAR·coverage·ASE·CNV는 기존 정책을 유지함. ASE 변이근거와 CNV에 필요한 WGS reads는 임시로 처리해야 하므로 `MFM_SEQUENCE_RNA_WGS_EVIDENCE=temporary`를 함께 지정해야 함. 지정하지 않으면 pair 실행 전에 실패함. 최종 WGS readset/depth는 저장하지 않음. WGS reads 자체를 전혀 사용하지 않는 RNA 정책은 별도 구현되지 않음.
- 읽기: `WindowReader(...).load(...)`는 서열·mask·RNA depth를 기존처럼 반환함. WGS depth만 `None`으로 반환하며 0 coverage로 해석하지 않아야 함.
- 검증·정리: 원본 FASTA와 모든 kept window를 독립 대조하고 배열·padding·N 주석·내용 SHA를 검증함. read QC와 구분한 `sequence_preflight.json`을 남김. 기존 최종 QC·source trace·소유권 확인을 통과한 뒤 FASTA/개인 VCF·임시 자료를 정리하고 서열·주석·압축 chain·QC·metadata·로그를 보존함. 실패 시 정리를 진행하지 않음.
- 자원: 기존 CPU/RAM 공유 풀을 사용함. 새 모드의 VCF 추출은 샘플별 CPU 위치를 분산하고 실제 worker 피크로 RAM을 반환하되 concat/검증용 최소4GiB를 유지함. CPU 수는 현재 grant를 넘겨 늘리지 않음. 서열 writer는 염색체별 multiprocessing이며 부모1GiB+worker당2GiB·64MiB 서열 배치로 RAM을 먼저 산정함. 제공 panel의 개인 VCF 추출은 최대22개 염색체 worker 기준20GiB, 개인 FASTA 단계는4GiB로 입장한 뒤 실제 worker 예산을 재확보함. 운영 풀의 실제 CPU/Slurm/cgroup·RAM 상한을 우선함.
- 경로: 모드가 다른 기존 완료 결과와 다른 `output_root`·`state_root`·`run_id`를 사용함. 기존 full-read 결과를 서열 전용 결과로 자동 변환하거나 덮어쓰지 않음. 기존 manifest의 용량 추정은 read 포함 기준일 수 있으며 실제 완료 용량은 저장 파일을 직접 집계함. RNA 출력은 DNA 서열 용량에 별도로 더해야 함.

기존 JSON을 복사해 `pipeline_env`에 다음 값을 지정한 뒤 내장 큐로 실행함. `--config`는 해당 JSON 경로이며 별도 운영 제어기는 필요하지 않음.

```json
"pipeline_env": {
  "MFM_SHORT_READ_OUTPUT_MODE": "sequence_only",
  "MFM_ZARR_BACKEND": "cpu",
  "MFM_MODEL_ONLY_CLEANUP": "1",
  "CLEAN_INTERMEDIATES_ON_SUCCESS": "1"
}
```

```bash
conda activate jj_mfm
bash run_pipeline.sh dry-run --selection manifests/리드별/medai/001/short_read.tsv --config /path/short_sequence.json --executor local
bash run_pipeline.sh run --selection manifests/리드별/medai/001/short_read.tsv --config /path/short_sequence.json --executor local
```

- 이 예시의 MedAI001 short104명은 WGS-only임. RNA pair가 포함된 다른 selection은 위의 임시 WGS 근거 설정도 필요함. `dry-run`은 선택·자원·예상 다운로드만 표시하며 실제 처리는 `run`에서 시작함.

## 재개 시 CPU·긴 read 메모리 보호 보완 (2026-10-28)

- CPU 자동 확장은 구실행 HAP 경계에서도 요청 변수 누락을 `auto`로 일관되게 전달함. 명시적인 고정 스레드 요청은 보존하고 호출 후 환경을 복원함. 이미 실행 중인 정렬기의 스레드를 실시간으로 바꾸는 기능은 아님.
- Compact read 필터는 `query_length`로 SEQ 부재를 판정해 불필요한 SEQ 문자열 캐시를 만들지 않음. 처리한 BAM 객체는 입력 순서를 유지하면서 즉시 해제함. read/window 필터·출력 순서·서열·depth 계약은 유지함.
- 길이 표본이 1,000bp 이상인 정렬은 worker 주소 공간 상한과 RAM 예약을 최소32GiB로 함께 설정함. 짧은 read의 하한8GiB와 입력 기반 더 큰 예측은 유지함. 전체1TB·예약900GB 안에서 동시 worker 수를 결정함. 이 하한은 극단적 국소 밀도의 절대 상한 보장이 아니며 실측 RSS 보호와 full/source QC를 함께 적용함.
- MedAI 큐는 실제 `EDQUOT`/`ENOSPC` 쓰기 오류 시 기존 작업 상태를 메모리에 보존하고 같은 큐에서30초 간격으로 재시도함. 여유 공간만으로 사전 중단하는 조건은 추가하지 않음. 검토된 실패 샘플은01/02/06/10/12단계에서 재개 가능하며, 부분 haplotag는 명시 요청한 경우 기존 백업 후 재생성 경로를 사용함.
- 이 서버 운영 정책은72코어 중 여유4코어를 제외한68코어 공유 풀임. 할당량과 실제 CPU 사용량은 구분함. RAM 상한·표준 full QC·독립 source trace·최종 QC 후 정리 정책을 유지함.
- VCF REF와 지정 reference 불일치는 그대로 실패 처리하며 임의 보정·검증 우회를 적용하지 않음.


## 완료 샘플 외부 원본·중간 정리 (2026-10-27)

- 모델전용 정리에서 샘플 내부 정리·정리 후 PASS 다음에 `completion_cleanup`을 실행함. 기존 실행도 이후 새로 시작하는 cleanup 서브프로세스부터 적용함. 별도 감시기 없이 기존 파이프라인에 연결함.
- 최종 QC·full compact·독립 source trace·정리 후 PASS와 필수 모델 입력을 확인한 뒤, 관리 저장소의 해당 샘플 전용 원본/인덱스·비공유 개인 VCF·외부 `_work`의 seq/tag/refok VCF·인덱스·part를 정리함. 공유 VCF/reference·다른 샘플·모델 입력·QC/metadata/log는 보존함.
- 열린 파일·symlink·계획 후 fingerprint 변경을 검사함. 외부 정리 실패는 별도 진행 보고서로 남겨 완료 판정을 막고, 모델 재생성 없이 외부 정리만 재시도 가능함.
- 기본 자동 정리는 소유권 등록 입력만 대상으로 함. 현재 완료25명에 대한 일회 정리에서는 이미 확인한 HG00189의 미등록 원본을 명시적으로 포함함. 향후 미등록 원본까지 일괄 자동 삭제하는 설정은 적용하지 않음.
- 2026-10-27 현재 완료25명의 잔여 데이터157개·67,565,196,278bytes를 제거함. 모델 검증 파일100개 해시 동일·삭제대상 잔존0 확인함. 작은 다운로드/출처 기록은 보존함. 감사자료는 `codex_archive/work/261027_완료외부정리/`에 있음.


## 현재 MedAI 경로·long pair 우선 배치 (2026-10-27)

- 실제 결과 루트: `/NHNHOME/korea/mfm_data/preprocessed_data/final/medai_001_261027`임. 다운로드 저장소는 그 아래 `raw/`임. `wgs/`·`wgs_rna_pair/`는 결과, `_work/`는 중간 작업 공간임. 완료 모델 입력은 최종 QC PASS 샘플의 `compact/`·`rna_personal/`만 선택함.
- 기존 `medai_policy_261020`과 `original/data`는 실행 중 작업·이전 절대경로 증거를 위한 호환 링크임. 실제 자료는 새 루트 한 곳에 저장하며 링크 경로를 중복 집계하지 않음. 현재 내부 N·blacklist 보존/주석·말단 N trim 정책을 유지함.
- MedAI 첫 논리 배치193명은 **HiFi long89명(pair10명 우선)+short104명**임. 002/003/004는각193명,005는192명의short임. 전체3212명·MedAI964명/TGIL2248명 기관배정은 유지함. TGIL 배치 변경·실제 실행은 없음.
- 현재 운영은 **long 전용·그중 pair 우선**이며 새 short는 시작하지 않음. 이미 실행 중인 short는 보존함. pair가 준비/자원 조건을 충족하면 먼저 투입하며, 다운로드·실패 대기 때문에 준비된 다른 long을 묶어 두지 않음.
- 기존 실행wave/state식별자는 별도로 보존함. 기존실행 `effective_config.assignment_snapshot`은 패키지 내부 `manifests/분담/운영호환_261023.tsv`의 SHA를 검증해 이전 배정에 대조함. 외부 selection URL을 신뢰하지 않고 동결 master 행을 사용함. 신규 `run_read_group.sh`는 이 호환설정을 제거하고 현행 배정표를 사용함.
- 재배정 도구: `tools/rebatch_long_first.py`임. 기본은 검토만 수행하고 `--write`로 목록을 갱신함. 인원·배치 용량 상한을 초과하면 실패함. 진행 중인 배정 변경 시 운영 state를 덮어쓰지 말고 검증된 snapshot과 별도로 관리함.
- CPU 실제 할당68코어·RAM1TB/예약900GB, 최종 full/source QC와 정리 후 검증·승인된 원본 정리 정책은 유지함.


## HiFi 단계 재배분·입력 병렬화 (2026-10-27 후속)

- **CPU 풀을 줄이지 않고, 각 단계가 실제 실행할 작업에 CPU와 RAM을 함께 예약함.** 이 서버의 실제 cgroup quota68코어·RAM 상한1TB·예약900GB 및 더 작은 Slurm 할당을 유지함. CPU 예약량과 실사용량을 구분함.
- 기존 단계 재배분이 빠졌던 **HiFi minimap2 경로에도 입력 집계 → target FASTA 준비 → 정렬 → index → 출력 QC → axis QC**를 연결함. 두 HAP 동시 실행은 각자의 현재 단계 요청을 합산함. 한 HAP가 QC로 넘어가도 다른 HAP의 정렬 예약은 보존함.
- minimap2의 mapping thread 외에 추가 I/O thread, view/fastq, sort main, unmapped awk를 포함해 **HAP당 보조CPU5개와 입력 decode thread**를 합산함. 신규 자동 설정은 sort 보조CPU 일부를 BAM 입력 압축 해제에 배분함. 실제 `samtools view -@ N`과 minimap2·sort의 합계가 grant를 넘지 않게 함. BWA는 기존 보조CPU3개를 유지함.
- 지원되는 public phased HiFi는 **현재 시작 단계 RAM으로 큐에 입장**함. 이미 정렬된 BAM의 Stage01은 실제 transcode 예약과 같은64GiB를 사용함. 개인 FASTA(원본09/legacy07)와 개인 재정렬은4GiB로 진입함. FASTA는 contig 작업당2CPU·4GiB와 부모4GiB를 확보한 뒤 실행하며, 재정렬도 도구 시작 전에 단계별 피크를 다시 확보함. 이전128GiB FASTA 계획이 환경에 남아 있어도 입장4GiB와 실제 worker 예약을 구분함. 전체 피크 모델 하한256GiB(WGS)·384GiB(WGS+RNA)는 단독 실행 가능성 검사와 계획 기록에 보존함. raw read 준비, 미검증 assembly/self-call 등은 해당 보수 모델을 유지함.

| 설정 | 기본값·동작 |
|---|---|
| `MFM_STAGE_ADMISSION` | `1`: 지원 HiFi의 시작 단계 예약을 사용함. `0`: 전체 피크 모델로 입장함. 후속 계산은 두 방식 모두 공유 풀에서 예약한 뒤 실행함. |
| `MFM_REALIGN_PHASE_RESOURCES` | `auto`: 지원 HiFi와 BWA에 단계 재배분을 적용함. `0`이면 종전 경로를 사용함. |
| `MFM_MINIMAP_DECODE_THREADS` | `auto`: HAP별 CPU 예산 안에서 최대4개를 선택함. `0` 또는 양수로 명시 가능함. 명시 요청이 전체 예산을 넘으면 실행 전에 차단함. |
| `MFM_REALIGN_INPUT_BACKEND` | `auto`: 단일 sample·32MiB 이상·유효 index·복수 CPU에서 구간 집계를 사용함. `serial`은 기존 전수 iterator, `indexed`는 index 기반 집계를 명시함. |
| `MFM_REALIGN_INPUT_CHUNK_BP` | `16777216`: 입력 집계의 구간 크기임. 양수로 변경 가능함. 모델 window와 무관하며 모델 입력 좌표·window 크기는 바꾸지 않음. |

- 입력 집계는 read의 **reference_start가 속한 구간에서만 한 번** 계수함. 경계를 가로지르는 read를 중복 집계하지 않음. 선택 염색체의 primary mapped, 모든 염색체의 primary unmapped, HP1/2/0·누락·이상값의 기존 규칙을 보존함. BAM index의 전체 record 통계와 범위 집계가 맞지 않으면 기존 전수 집계로 복귀함. CRAI의 unmapped 통계 한계는 독립 native 전수 count로 확인함.
- 집계 worker당 RAM2GiB와 부모4GiB를 예측하고 CPU·RAM·작업별 주소 공간 한도를 함께 지킴. 집계의 임대를 임의로 축소해 재배분이 꺼진 후속 정렬이 예약 없이 실행되는 경우를 막음. NM/axis QC의 decoder도 해당 하위 단계 grant에 맞춤.
- 다운로드의 전체 BAM/CRAM CRC 검사 결과를 최종 경로의 `.download_records/<key>.integrity.json`에 연결함. 같은 inode·크기·mtime의 통제된 `.part` rename만 인정함. worker는 최종 fingerprint와 동일 scanner SHA를 확인해 공유 결과를 재사용하며 자체 검증 보고서에도 남김. 다운로드 완료 게시 전에2초 안정화와 fingerprint를 확인함. worker에서는 기존 파일-age 보호를 대기로 대체하지 않고, 파일 변경·손상 시 원래 전수 재검증을 유지함. **검사를 완화하는 방식이 아니라 같은 파일의 중복 전수 읽기를 제거하는 방식임.**
- 실물 소형 HiFi BAM에서 종전/신규 minimap2의 read·tag·QC가 같고, BAM/CRAM의 배치·미배치 unmapped와 HP 분류가 같음을 확인함. 실행 중인 mapper의 thread는 소급 변경하지 않으며 새 설정은 새 단계 프로세스부터 적용함. NA18959는 HAP1/2 정렬 QC를 통과해 저장 단계로 진행했고, NA20850의 기존 정렬은 보존함. 준비된 NA21093도 추가 투입됨. 04:21–04:22의30초 표본은 실제49.48코어/풀68코어(예약67코어)임. 서로 다른 단계의 이전25.1코어 표본과 비교한 전체 속도 배수나 상시68코어 실사용을 뜻하지 않음. 원본125건·전달용120건 회귀, 공유 무결성 캐시 각7건, 운영큐61건 검사를 통과함.

## CPU 유휴 구간 완화·HiFi 예약 개선 (2026-10-27 추가)

| 개선 항목 | 현재 동작 |
|---|---|
| 샘플 사이 대기 | `queue_policy.wait_for_drain=false`이면 기존 실행을 보존하면서 준비된 long을 빈 자원에 즉시 투입함. `long_read_only`에서는 short 계산·다운로드로 자동 복귀하지 않음. |
| 큐 RAM 판단 | 실행 샘플의 생존 stage lease를 사용함. 일괄100GiB/sample 하한·지난 index 피크의 중복 예약을 제거함. lease가 없으면 보수적 전체 샘플 모델로 복귀함. |
| HiFi 전체 피크 모델 | public phased VCF + `vcf_consensus`의 준비·VCF/consensus·개인 정렬/sort/QC·RNA 피크를 reference 길이·최대 thread·명시 버퍼로 예측함. 입력 파일은 streaming 처리하며 크기4배 적재를 가정하지 않음. |
| 기존 순차 HAP | 명시 요청한 작업만 HAP1 정렬·index·QC가 모두 끝난 뒤 HAP2 실행 직전에 재배정함. 현재 BWA는 유지하며 완료 HAP를 다시 처리하지 않음. |
| 정렬 후 예약 반환 | 순차 HAP의 alignment pipeline 종료 뒤 index·출력 QC·axis QC의 CPU/RAM을 각각 다시 예약함. 두 HAP 동시 계산 중인 파이프라인의 보조CPU+3/HAP는 유지함. |
| 계산 내부 병렬화 | FASTA는 HAP×chrom 독립 작업으로 분할함. full QC의 전체 read CIGAR 검증은 bounded NumPy chunk로 수행하며 남는 CPU/RAM만 추가 lease로 빌림. window depth 재검증·해시·독립 source trace는 유지함. |

- **전체 RAM 상한1TB·예약900GB 및 더 작은 Slurm/cgroup 할당을 우선함.** CPU 설정72인 이 서버의 실제 cgroup quota는68코어임. 코드가 동시 샘플 수 때문에 풀을 축소한 값이 아님. CPU 예약량과 실제 계산 사용량을 구분함.
- `MFM_LONGREAD_MEMORY_MODEL=stages_v1`이 기본임. GRCh38/decoy reference 상한40억bp, 최대CPU/thread/sort buffer와10%여유를 적용함. 전체 피크 모델의 하한은 일반 WGS256GiB, WGS+RNA384GiB임. 현재 시작 예약은 위의 단계 입장 정책을 적용함. 공통 부모·다운로드 예산은 큐에서 별도 합산함. 실행 중 compact는 실제 read/CIGAR 수로 RAM을 다시 예약함.
- 같은 풀에 과거 할당기가 남아 있어도 새 도구 실행 직전에 요청 RAM과 실제 lease를 다시 대조함. 다르면 남은 CPU/RAM 안에서 재예약하고, 부족하면1CPU·4GiB 대기로 반환함. 이미 시작한 자식의 예약은 이 보호로 바꾸지 않음.
- `MFM_LONGREAD_MEMORY_MODEL=legacy`로 종전 입력크기 기반 예약을 선택 가능함. 자체 변이호출·assembly·비표준 profile 등 미검증 경로는 legacy를 유지함. 실제 메모리 피크 측정값으로 해석하지 않음. `ram_admission_plan.json`과 sample `metadata/longread_memory_plan.json`에 계획을 보존함.
- 1TB 조건에서 새 모델을 충족하는 MedAI long89명·TGIL long205명은 개별 시작 예산 이내임. 실제 동시 실행은 현재 lease·RSS·CPU/RAM 여유에 따라 대기 가능함. `manifests/리드별/RAM재검토.tsv`는 과거19명 추정값을 보존하고 새 모델 열을 추가한 비교표임. TGIL 목록 검증은 TGIL 서버 실실행을 의미하지 않음.
- `MFM_FASTA_CONTIG_JOBS=auto`가 기본임. 단일 sample 호출에서 HAP×chrom별 최대44태스크를 CPU2개·RAM4GiB/task 범위로 실행함. 숫자1은 직렬임. 여러 sample 직접 호출은 기존 sample 병렬을 유지하고 내부 중첩 병렬을 막음. FASTA·chain·QC 결과 순서를 유지함.
- `MFM_QC_CIGAR_BACKEND=numpy`와 `MFM_QC_CHUNK_THREADS=auto`가 기본임. 최대262,144레코드·8,388,608 CIGAR op 단위로 자름. 너무 큰 범위는 재분할하며 단일 거대 read는 기존 decode로 검증함. 추가 thread당2GiB를 예약하고 주소 공간 여유도 확인함. Blosc가 이미 여러 thread이면 추가 QC thread를 사용하지 않음. backend `legacy` 또는 thread1로 비교 가능함.
- 기존 작업의 선택적 경계 요청은 공유 pool의 `hap_boundary_requests.json`을 사용함. manager PID/시작시간·소유 경로·HAP1 QC/index·생성 shell 형식·FAI를 검증함. HAP2 shell은 백업하고 숫자 thread만 원자 변경하며 `HAP2.boundary_upgrade.json`에 실제 설정을 보존함. 요청하지 않은 작업은 기존 동작을 유지함.
- public VCF index가 다운로드 목록에 없는 경우 사전 검증 직전에 누락 CSI만 생성함. 경로별 잠금·임시 index 검증·원자 교체를 사용하며 원본 VCF와 기존 index는 보존함. 실제 stage grant 이내 최대4CPU를 사용함. Stage01 완료 후 index 누락으로 실패한 작업은 검토 후 Stage02부터 재개 가능함.
- 작은 실제 VCF에서 순차·병렬 FASTA/chain의 바이트와 QC가 같음을 확인함. 실제 BWA·BAM/CRAM·Zarr 손상 검증·lease 반환/상한 검사를 수행함. 전체 HiFi 전장 처리시간 감소율과 CPU68코어 지속 사용은 아직 검증하지 않음.

## 정렬 단계 CPU·RAM 재배정 (2026-10-27)

- Illumina BWA-MEM2 정렬의 RAM 예약은 `reference_v2`를 기본으로 사용함. reference 길이×16bytes, HAP별24GiB의 기본/QC 버퍼, BWA thread당256MiB, HAP별4GiB 파이프라인 버퍼, sort의 `-m × (보조thread+1)`을 합산하고10%여유·최소32GiB·GiB단위올림을 적용함. **전체1TB상한·예약900GB를 유지하며, 실측RSS만으로 예약을 낮추지 않음.**
- index 생성은 별도40bytes/reference-base 예약을 유지함. 입력 준비·HiFi 전체샘플 예약과 정렬 RAM을 같은 지표로 해석하지 않음. 정렬모델만으로는 HiFi 보류를 해제하지 않으며, 위의 별도 단계별 모델 검증을 적용함.
- `MFM_REALIGN_MEMORY_MODEL=legacy192`로 기존192GiB하한 모델을 선택 가능함. 명시thread/hap요청은 유지함. 모델의 근거·구성요소는 `effective_alignment.memory_model`과 resource phase 기록에 남김.
- 순차HAP는 앞HAP 정렬·index·QC 완료 뒤 새HAP 시작 직전에 CPU/RAM을 다시 계산함. 병렬HAP 중 하나가 완료되면 그 몫의 예약을 반환하고, QC도 남은 실제HAP 수를 사용함. 실행중BWA의 `-t`는 변경하지 않음.
- 일반 경계 재배정은 새 전처리Python부터 적용함. 이미 실행 중인 BWA는 유지함. 명시 요청한 기존 순차 작업은 위의 resource_runner HAP2 경계 연결을 사용하며, 아래 도구는 검증된RAM예약만 줄임.

```bash
# 기본은계획만출력함. 실제검토후 같은명령에 --apply를추가하면예약만감소함.
python tools/rebalance_realign.py \
  --pool-dir /실제/공유CPU풀 \
  --cohort-root /실제/출력/wgs/1KGP_ILLUMINA \
  --sample SAMPLE_ID --report /실제/기록/예약계획.json
```

- 보정은 실행중BWA·manager/PID시작시간·실제thread·개인FAI·sort버퍼·RSS를 확인하고 같은pool.lock 안에서만 수행함. 새모델과실측RSS두배 중 큰 값을 예약하며, 감소할 때만 허용함. index/대기/이미끝난 작업, lease/PID/reference 변경, RSS여유 부족은 거부함. CPU/출력파일/주소공간제한을 수정하지 않음.
- 최종full QC·독립source trace·성공후정리는 유지함. 작은실제BWA에서 HAP경계7→17thread재배정과 read/QC동등성을 검증함. 전장속도개선배수는 단계변화와구분해 별도측정해야함.

## Long/short read별 실행 목록 (2026-10-26)

`manifests/리드별/`에 두 기관의 전체 종류별 목록과 기존16개 배치별32개 TSV를 포함함. **MedAI long89·short875명, TGIL long205·short2043명**임. long은 HiFi, short는 Illumina이며 파일 확장자로 구분하지 않음. 기존 기관/배치·원본 manifest 칼럼과 전처리 정책을 유지함.

```bash
# 실행 서버에 맞춰 기존 예제의 data_root/output_root/state_root/자산 경로를 지정함.
# jj_mfm 활성화 후 이 패키지 폴더에서 실행함.
bash ./run_read_group.sh dry-run --lab tgil --read-type long_read --batch all --config handoff_config.example.json
bash ./run_read_group.sh run --lab tgil --read-type long_read --batch all --config handoff_config.example.json

# 단일 배치 또는 short read도 같은 방식으로 선택함.
bash ./run_read_group.sh dry-run --lab tgil --read-type short_read --batch 001 --config handoff_config.example.json
```

- `--lab medai|tgil`, `--read-type long_read|short_read`, `--batch 001|all`을 사용함. `all`은 기존 배치별로 순차 실행하며 한 실행에서 기관/배치를 섞지 않음. local 실행용 보조 진입점임.
- 내부에서 기존 `run_pipeline.sh run --selection ... --executor local`을 호출함. 원본 다운로드→WGS/RNA→full QC/source trace→성공 후 소유 원본·중간 파일 정리 정책을 기존 config대로 적용함. long read도 기본 vcf_consensus를 유지하며 assembly 직접입력으로 바꾸지 않음.
- 전체 목록은 `manifests/리드별/tgil_long_read.tsv`, `tgil_short_read.tsv`, `medai_long_read.tsv`, `medai_short_read.tsv`임. 실제 단일 배치 선택은 `manifests/리드별/tgil/001/long_read.tsv` 등의 경로임.
- 출력은 `<output_root>/<lab>_<batch>/...`를 유지함. 선택 지문 충돌을 막기 위해 상태만 `<state_root>/read_groups/<read_type>/<lab>_<batch>/...`로 구분함. 같은 설정의 long/short는 동일 CPU/RAM 공유 풀을 사용함. 서로 다른 제어기에서 같은 출력 샘플을 중복 실행하지 않음.
- 이미 서버 전용 운영 큐가 실행 중인 MedAI001에는 보조 진입점을 추가 실행하지 않음. 기존 큐의 고정 선택/task_index를 그대로 유지하고 별도 우선순위 정책으로 순서만 전환함.
- 보조 실행은 실제 RAM 상한의90% 예산과 공통16GiB·다운로드12GB를 포함해 단독 투입 불가 항목을 먼저 분리함. 제외 사유는 batch state의 `launch/제외사유.json`, 전체 요약은 `read_groups/<lab>_<read_type>_summary.json`에 기록함. 보류가 남으면 `complete_with_pending`·종료1로 보고하며 완료로 세지 않음. 예상 예약량이며 실제 RAM 피크와 구분함.
- 종전 legacy 산식의1TB 예비 분류는 MedAI long81명 가능/8명 재검토, TGIL long194명 가능/11명 재검토였음. 현재는 위 stages_v1 모델로 각각89명·205명의 개별 예산을 다시 검토함. 더 작은 Slurm/cgroup/사용자 예산에서는 대상이 더 줄어들 수 있음. `manifests/리드별/RAM재검토.tsv`는1TB 기준 정적 참고표이며 실제 실행 시 재계산함.
- 기본 큐는 준비완료 우선임. 신규 혼합 배치에서 long 우선을 원하면 config에 `"queue_policy": {"mode": "long_read_first"}`를 지정함. 해당 배치에서 처리 가능한 long이 완료/실패보류되기 전까지 short 신규 계산/다운로드를 보류함. 준비된 short 파일은 보존함.
- short 신규 계산/다운로드를 계속 중지하려면 `"queue_policy": {"mode": "long_read_only"}`를 사용함. long이 모두 완료되거나 보류돼도 short로 자동 복귀하지 않음. 기존 실행 샘플을 끝까지 보존할 때 해당 ID를 `drain_samples`에 명시함. 중지한 short는 실패 장부에 추가하지 않고 `priority_policy.paused_short_samples`로 구분함.
- 운영 중 정책은 `<queue 상태경로>/우선처리정책.json`으로 지정 가능함. `mode`, 기존 작업을 보존할 `drain_samples`, 전체 완료를 기다릴지 정할 `wait_for_drain`(생략 시 true), 사유를 명시한 `blocked_samples`를 사용함. 기존 실행 작업은 취소하지 않으며 long 선다운로드는 drain 중에도 가능함. 자원 상한 때문에 임의로CPU풀을 줄이지 않음.
- 동시 샘플 수는 같은 정책의 `"max_active_samples": 4`로 실행 중 변경 가능함. 양의 정수는 신규 계산 입장 상한이며 `null`은 고정 상한 해제임. 생략하면 시작 시 `--local-jobs`/자동 설정을 사용함. 상한을 초과해 이미 진행 중인 샘플은 중단하지 않고 완료까지 보존하며, 활성 수가 상한보다 작아질 때만 다음 샘플을 투입함. CPU 풀·RAM 예산·선행 다운로드 수는 변경하지 않음.
- 기존 `model_target_bytes`는 유지함. 목표 도달 시 `target_reached`로 멈추며 모든 샘플이 끝났다는 의미와 구분함. 실측 전장 처리시간·TGIL 서버 실실행은 이 목록 검증에 포함하지 않음.

## 현재 권장 실행과 자동 큐 (2026-10-29)

원본과 전달용 모두 `run_pipeline.sh run`의 **local 실행은 공용 연속 큐**를 사용함. 실행 코드는 날짜·서버에 종속된 운영 Python을 import하지 않음. 기존 실행 장부는 최초 한 번 호환 정보로 연결 가능함. 원본의 기존 `sample`, `wgs`, 단계별 다운로드 명령도 유지함.

```bash
# 예제 config를 복사해 data_root/output_root/state_root를 실행 서버에 맞게 지정함.
bash ./run_pipeline.sh dry-run --lab medai --batch 001 --config handoff_config.example.json
bash ./run_pipeline.sh run --lab medai --batch 001 --config handoff_config.example.json --executor local
```

| 항목 | 현재 동작 |
|---|---|
| 계산 투입 | 다운로드가 끝난 샘플부터 즉시 투입함. 묶음 전체의 완료를 기다리지 않음. |
| 상태 확인 | 1초 주기임. 준비·검증·자원 대기까지 1초 안에 끝난다는 의미는 아님. |
| `local_jobs=auto` | 확정한 자동 동시수 안에서 RAM·cgroup·RSS로 신규 입장을 결정함. 숫자 명시 시 해당 동시 수 상한을 유지함. |
| `download_jobs=auto` | 최소3개 또는 현재 계산 수에 맞춰 선다운로드함. 대기/다운로드 버퍼는 기본 목표의2배임. 숫자 명시 시 해당 수를 사용함. |
| `prefetch_samples` | config에서 버퍼 상한을 따로 지정 가능함. 미지정 시 다운로드 목표의2배임. |
| 실패 | 다운로드는 샘플별 최대3회 재시도함. 지속 실패·계산 실패는 해당 샘플만 보류하고 다른 샘플은 계속 진행함. |
| CPU | 실제 affinity/cgroup/Slurm 할당을 우선함. 동시 샘플 수 때문에 전체 CPU 풀을 줄이지 않음. |
| RAM | 운영 상한1TB·예약90%이며 더 작은 할당을 우선함. 불명확한 단계는 보수적 예약을 유지함. |
| 디스크 | 공용 local 큐는 가용량으로 신규 투입을 차단하지 않음. 제어 큐 쓰기 오류 ENOSPC/EDQUOT는30초 간격 재시도함. 개별 worker 실패는 해당 샘플만 보류함. |
| 완료·정리 | 기존 full QC·독립 source trace·정리 후 검증 계약을 유지함. 완료본은 원본이 없어도 재다운로드하지 않음. |
| 모델 용량 목표 | 완료 모델 입력이 `model_target_bytes`에 도달하면 신규 투입을 멈추고 시작한 샘플은 마무리함. 동시 실행분 때문에 목표를 초과할 수 있음. |

현재 예제 config는 모델 입력 전용 정리(`MFM_MODEL_ONLY_CLEANUP=1`, `CLEAN_INTERMEDIATES_ON_SUCCESS=1`)와 소유 원본 정리를 명시함. 성공 검증 전에는 정리하지 않으며 공용 reference는 보존함. 별도 사용자 config의 값은 보존함. 현재 서버의 기본 MedAI config는 모델 전용 정리·CPU풀68·long-only/pair 우선을 명시함.

명시한 `WINDOW_STORAGE_POLICY`와 기존 완료본의 정책이 다르면 그 결과를 건너뛰거나 덮어쓰지 않고 별도 output/state 경로를 요청함. 기존 선택 fingerprint와 목록 순서가 달라도 상태 저장 전에 거부함.

최초 공통 reference 준비는 계산 전에 수행함. 데이터별 입력 검증·CPU/RAM 부족·원격 응답 지연이 있으면 대기할 수 있음. `run --executor slurm`은 기존 prepare→array→finish 제출 방식이며, local의 독립 선다운로드 큐가 Slurm array에도 적용된다고 보장하지 않음.

상태와 재시도·대기 이유는 `<배치 state_root>/queue/status.json`, `events.jsonl`, `재시도상태.json`에 저장함. `queue/control.lock`으로 중복 제어기를 차단하고, 재시작하면 기존 worker/download의 잠금을 채택함. 종료 요청을 받은 제어기는 새 투입을 중지하고 자신이 시작한 worker가 끝날 때까지 대기함. 실패한 계산은 아래 `run_pipeline.sh resume` 명령으로 재개 요청함. 필요한 JSON은 코드가 자동 생성함. 준비 파일과 원본이 남아 있는 검토된 단계에만 적용함.


### 단일 진입점과 CPU 자동 계산 (2026-10-29)

- 기본 설정: `handoff_config.json`을 사용함. 파일이 없으면 `handoff_config.example.json`을 사용함. 외부 전달 서버의 경로/대상을 처음 한 번 설정한 뒤 같은 명령으로 재개 가능함. `--config FILE`은 기본 설정을 바꿀 때만 필요함.
- 기본 대상: config에 `default_selection`을 지정하면 `bash run_pipeline.sh` 또는 `bash run_pipeline.sh run`만으로 해당 목록을 실행함. 미지정 서버는 `--lab/--batch` 또는 `--selection`을 명시해야 함. 무인자로 전체 3,212명을 실행하지 않음.
- 자동 실행: 내장 `handoff/local_queue.py` → `ready_queue.py`가 선다운로드·완료 입력 우선 투입·자원 검사·실패 격리·재시작 인수를 수행함. 별도 `연속제어.py` 등 운영 Python을 사용자가 실행하거나 수정할 필요 없음. 실제 샘플은 기존 `run_sample.sh`로 실행함.
- CPU 풀: `C = min(affinity, cpu_pool_limit, cgroup CPU 할당, Slurm 실제 할당)`임. 미지정/`auto`인 상한은 제한 요소에서 제외함. `cpu_utilization=maximum`이면 동시 샘플 수와 관계없이 이 풀을 공유함. `bounded` 로컬 모드는 추가로 `sample_cpus × local_jobs` 이하로 제한함.
- 초기 자동값: `sample_cpus = max(1, C // min(J, 4))`임. J는 profile·가용 RAM(샘플당 128GiB 계획)·CPU 조건으로 계산한 초기 `local_jobs`임. C=68일 때 J=1/2/3/4 이상이면 각각68/34/22/17임. local 큐의 `local_jobs=auto`는 확정한 J를 신규 입장 상한으로 사용하며 실제 단계 RAM·cgroup·RSS 검사를 함께 수행함. 이미 진행 중인 계산은 보존함.
- 실제 배분: 단계별 최소/최대 CPU·예상 RAM을 요청함 → 생존 lease를 제외한 빈 CPU를 준비된 요청끼리 나눔 → 단계 종료 시 반환함. 긴 계산은 최대 공유 풀까지 요청 가능함. 이미 실행한 정렬기의 thread 수는 강제로 바꾸지 않음. 지원되는 내부 worker는 경계에서 남는 CPU/RAM을 재요청함. 직렬 단계와 I/O 구간은 100% CPU 사용을 보장하지 않음.
- RAM: 운영 상한1TB·예약900GB 및 더 작은 실제 할당을 우선함. 선다운로드는 전송당4GB를 검사·예약함. 디스크 잔량으로 미리 중지하지 않음.
- 저장 오류: 제어 큐의 실제 `ENOSPC`/`EDQUOT` 오류는 기존 worker 객체·메모리 장부를 보존한 상태로30초 후 재시도함. 개별 계산 자체가 실패한 경우에는 자동으로 임의 단계를 선택하지 않고 검토된 단계 재개 요청을 사용함.

```bash
# 기본 설정/대상이 등록된 서버: 이 명령만 실행함.
bash run_pipeline.sh
# 설정 결과 확인: 다운로드/전처리 없이 확인함.
bash run_pipeline.sh dry-run
# 전체 CPU 공유 풀 상한을 명시함.
bash run_pipeline.sh run --cpu-pool-limit 68
# 실패 샘플의 검토된 단계부터 재개 요청함. JSON은 자동 생성함.
bash run_pipeline.sh resume --sample SAMPLE_ID --from-stage 01
```

- 재개 단계: 원본01/02/06/10/12, 전달용01/02/06/08/10 허용함. 준비 파일과 원본이 필요한 상태인지 먼저 확인함. 활성 큐가 있으면 다음 확인 주기에 적용함. 큐가 중지된 상태이면 같은 `run` 명령으로 이어감. 실행 중인 샘플은 요청을 거부함. `--overwrite-prepared`·`--rebuild-haplotag-partial`은 검토 후 명시할 수 있음.
- 장부: `effective_config.json`, `selected_samples.tsv`, `samples/`, `queue/status.json`, `queue/events.jsonl`, `queue/재시도상태.json`, `queue/작업재개요청.json`을 코드가 자동 생성함. 기본적으로 사용자가 직접 작성하는 파일이 아님.
- 기존 운영 이관: `run_pipeline.sh import-state --from-dir 이전장부폴더`를 한 번 사용함. 이전 제어기의 신규 투입을 중지한 뒤 실행함. 선택/출력 경로를 검증하고 과거 state·task index·완료/실패/재개 장부를 연결함. 외부 운영 Python을 실행·import하지 않음. 활성 worker의 과거 state·공유 pool 경로는 계산 종료까지 보존해야 함. 새 서버의 신규 실행에는 이관 명령이 필요 없음.

### 선택 GPU 백엔드 (2026-10-29)

- 기본값: 기존 CPU 처리 유지함. WGS·RNA·Coverage·QC·Zarr 패킹을 각각 선택 가능함. 동일 `run_pipeline.sh` 진입점과 기존 다운로드→전처리→full QC/source trace→성공 후 정리 흐름을 사용함.
- 현재 검증: B200에서 정수 depth·RNA deletion/splice/duplicate 필터·4bit bytes/offset·CIGAR 오류·생성 Zarr의 full semantic QC·독립 CPU source trace 및 변조 탐지 검증함. 실제 Parabricks 정렬의 전장 CPU/GPU 일치성·처리시간은 별도 실험 대상임.

| 항목 | CPU 선택 | GPU 선택 | 실제 GPU 계산 및 제한 |
|---|---|---|---|
| WGS alignment | `--wgs-backend cpu` | `--wgs-backend parabricks` | HiFi/ONT는 GPU minimap2, Illumina는 GPU BWA-MEM 사용함. HAP1/2별 정렬·기존 routing·OD duplicate 복원·최종 CRAM/QC 유지함. CPU BWA-MEM2 index와 GPU용 BWA-MEM index는 분리함. |
| RNA alignment | `--rna-backend cpu` | `--rna-backend parabricks` | 별도 호환 STAR2.7.2a로 개인 index 생성 후 GPU STAR 사용함. 기존 CPU STAR2.7.11b와 구분함. |
| Coverage | `--coverage-backend cpu` | `--coverage-backend cupy` | CPU에서 BAM/CIGAR를 해석하고 GPU에서 구간별 scatter-add·prefix sum을 수행함. WGS는 M/=/X depth, RNA는 deletion 포함·splice N 제외·MAPQ255·duplicate 포함·mate 겹침 두 번 계수함. |
| QC | `--qc-backend cpu` | `--qc-backend cupy` | CIGAR 유효성·query/reference 소비 길이·정수 depth 통계를 GPU에서 계산함. 독립 source trace와 semantic depth/sequence 대조는 CPU 구현을 보존함. |
| Zarr 준비 | `--zarr-backend cpu` | `--zarr-backend cupy` | 4bit LUT 변환·유효염기 검사·패킹을 GPU에서 수행함. 최종 저장은 기존 Zarr2·Blosc/Zstd5·uint32 depth 유지함. 압축·파일 쓰기를 전부 GPU로 수행한다는 의미는 아님. |

### 앞단 GPU IO 선택과 검증 (2026-10-29)

- `--io-backend cpu|nvcomp` 선택 추가함. **현재 기본/운영은 CPU임.** GPU IO는 구현·동등성 검증했으나 이 서버 실측에서 CPU보다 느려 자동 적용하지 않음.
- GPU IO 범위: CRAM3 block IEEE CRC32 배치 계산과 canonical BAM의 BGZF 압축/쓰기임. 원본 BAM BGZF 손상 검출·MD5·CRAM 해독·samtools comparator 정렬·index·VCF REF 확인은 CPU 유지함. `--qc-backend cupy`와 별도 설정임.
- 정확성: GPU 압축 모든 chunk를 CPU로 독립 해독하여 원본 바이트·CRC와 대조함. 실제1.901GB BAM 전체 SHA256 동일함. Stage01 CPU/GPU 변환의 record/tag/순서와 CRAM3.0/3.1 손상 검출 offset 일치함.
- 성능: 같은4CPU 예산의 NA18959 실제 read 압축/쓰기 CPU7.220초/GPU8.058초임. 출력784.4MB/906.0MB임. HG00275 원본14.807GB 전체 물리 CRC CPU13.471초/GPU15.580초임. 전체 전처리 가속률이 아닌 IO 구간의 실측임.
- 의존성: 선택 requirements-gpu.txt의 `nvidia-libnvcomp-cu13==6.0.0.12` C API 사용함. 설치는 기존 `00_environment.sh install-gpu` 사용함. 이 IO 구현은 CUDA13 nvCOMP6 ABI를 검증했으며 다른 변형은 별도 검증 필요함.
- 자원: GPU IO 자식이 공용 GPU 풀의4GiB 예약을 사용하고 종료 후 반환함. bounded 미리읽기2배치/독립 검증2배치 사용함. canonical GPU IO는6CPU 이상에서 사용하고 더 작은 grant는 CPU로 수행함. CRAM 검증 자식은 별도 부모 stage 예약이 없으면 IO CPU2/RAM4GiB를 공유 풀에서 예약함.
- 저장/기록: BAM 압축 바이트·크기·@PG는 달라질 수 있으나 read/SEQ/QUAL/CIGAR/FLAG/태그/정렬 순서를 보존함. 실제 GPU IO는 Stage01 QC의 source_io_execution과 stage 로그에 기록함. 최종 모델 형식/윈도우/필터 정책·CPU 독립 source trace는 유지함.

**RNA GPU의 정책 보호:** 기본 strict에서는 현재 GTEx profile의 GPU 실행을 차단함. 차이를 평가하는 실험은 `--rna-gpu-policy experimental`을 명시해야 함.

**RNA GPU의 실제 정책 차이:** 현재 Parabricks4.8은 `outFilterType=BySJout`을 설정하는 인터페이스를 제공하지 않아 GPU STAR 기본 `Normal`을 사용함. STAR2.7.2a와 현재 CPU2.7.11b의 차이도 존재함. GPU 결과에는 별도 profile 이름·실효 인자·비동등성·버전·이미지·GPU UUID를 기록하고 재사용 지문을 구분함. 기존 GTEx CPU 정책과 동일한 결과로 취급하지 않음. 실험 모드를 명시한 비교 실험에서 먼저 평가할 대상임. 기존 CPU 정책은 유지함. 또한 CPU의 `limitSjdbInsertNsj=1200000`은 GPU에서 기본1000000임. `max-collapsed-junctions`는 다른 설정인 `limitOutSJcollapsed`이므로 이를 삽입 한도의 대체 인자로 사용하지 않음. 실제 index의 `sjdbOverhang`을 GPU 정렬에도 명시함.

```bash
# 설정/명령 검토: 실제 다운로드나 GPU 계산을 시작하지 않음.
bash run_pipeline.sh dry-run \
  --selection manifests/리드별/medai/001/long_read.tsv \
  --wgs-backend parabricks --rna-backend cpu \
  --coverage-backend cupy --qc-backend cupy --zarr-backend cupy \
  --gpu-devices auto --gpu-runtime auto

# GPU 의존성 및 호환 CPU 도구 준비함. 기본 core STAR/BWA-MEM2를 교체하지 않음.
bash 00_environment.sh install-gpu
bash 00_environment.sh install-gpu-tools
bash 00_environment.sh verify-gpu
```

- 정렬 실행환경: NVIDIA Parabricks4.8.0-1 이미지와 Docker/Apptainer, 또는 준비된 native `pbrun`이 필요함. `--gpu-runtime auto|docker|apptainer|native`, `--parabricks-image IMAGE`, `--gpu-rna-star FILE`로 지정 가능함. CuPy 설치만으로 Parabricks가 설치되는 것은 아님.
- 실행환경 갱신(2026-10-29): 공식 NGC4.8.0-1 application 레이어를 SHA256 검증 후 사용자 경로에 준비하여 native 실제 WGS 실행 검증함. `MFM_PBRUN`에 준비된 실행 파일을 지정하고 `--gpu-runtime native`로 선택 가능함. 기본 core 도구와 시스템 설정은 교체하지 않음. 현재 MedAI 서버의 HG00275·NA18959·HG00126 전체 샘플 대조는 GPU0·1·2 공용 풀과 별도 output/state에서 진행 중임. 다른 서버의 native 호환성이나 전장 가속률까지 검증된 것은 아님. `install-gpu-tools`는 STAR/BWA 호환 도구 설치이며 Parabricks runtime 자동 설치 명령은 아님.
- 실제 실행: 위 `dry-run`을 `run`으로 바꿔 실행함. 이미 동일 큐가 실행 중이면 중복 제어를 차단함. CPU/GPU 결과를 비교할 때에는 별도 config의 state/output 경로를 사용하여 기존 완료본과 구분함. 기존 완료본을 자동 덮어쓰지 않음.
- 자원: CPU68·host RAM1TB/예약900GB 및 더 작은 실제 할당을 우선함. `CUDA_VISIBLE_DEVICES`·`NVIDIA_VISIBLE_DEVICES`와 명시한 `--gpu-devices` 범위 안에서 선택함. 외부 사용 흔적이 있는 GPU는 신규 선택에서 제외함. 정렬은 GPU 한 장을 독점 예약하고, CuPy worker는 기본4GiB씩 VRAM을 예약해 공유함. 공유 계산을 한 장에 우선 모아 다른 GPU를 정렬용으로 남김. GPU 예약은 이 파이프라인끼리의 협력적 잠금이며 외부 작업의 GPU 사용을 강제로 막는 기능은 아님.

- GPU 공용 풀 실행: 모든 샘플에 같은 `--gpu-devices 0,1,2`와 `gpu_pool_dir`를 지정하면 준비된 샘플/HAP 작업이 빈 GPU를 동적으로 예약함. 샘플별 GPU 번호와 앞 샘플 종료 대기는 필요하지 않음. `pipeline_env`의 `CUDA_VISIBLE_DEVICES`도 동일 목록으로 두고 `MFM_GPU_HAP_JOBS=2`, `WGS_HAP_JOBS=auto`로 설정하면 HAP 동시 수는 CPU/RAM 예산으로 결정함. 정렬은 HAP당1장을 독점 예약하고 완료 후 풀로 반환함.
- 단계별 변경 반영: 내장 `cpu_pool.py`가 같은 run/input/output의 `effective_config.json`에서 GPU 목록·공유 예약 경로·HAP 설정만 다음 단계 시작에 다시 읽음. 이미 실행 중인 도구는 유지함. backend/입력/출력 설정의 실시간 변경은 지원하지 않음. Slurm CUDA 공개 범위와 `NVIDIA_VISIBLE_DEVICES`는 확장하지 않음.

- HAP 동시성: GPU WGS는 기본 한 GPU에서 HAP1/2 순차 실행함. 복수 GPU가 허용된 환경에서는 `MFM_GPU_HAP_JOBS=2`로 선택 가능함. GPU host 계산은 기본 HAP당16 CPU 범위에서 계획하고, 정렬 도구 host RSS도 별도 감시함. GPU의 큰 가상 주소 공간은 실제 RAM과 구분하며, CuPy 초기화 후 AS 한도를 재산정하고 공유 RSS 보호를 유지함.
- Slurm: GPU 모드가 하나라도 선택되면 sample array에 기본 `--gpus-per-task=1`을 요청함. 준비·요약 작업은 GPU를 요청하지 않음. 필요하면 config의 `slurm.gpus_per_task`와 GPU partition을 설정함. 실제 Slurm GPU 제출은 현재 서버에서 수행한 검증에 포함하지 않음.
- 기록: `effective_config.json`의 acceleration/pipeline_env, WGS `logs/*.gpu.json`, RNA 최종 `rna_realign.json`의 gpu_execution/재사용 지문, Zarr `compute_execution`에 기록함. 중간 GPU FASTQ/BAM은 작업용 임시 경로를 사용하고, 최종 모델 입력의 dtype·padding·좌표·필터 기준은 유지함.
- 환경 선언: `requirements-gpu.txt`는 CUDA13.x용 CuPy14.2.0을 고정함. CUDA12.x 환경은 CuPy를 중복 설치하지 않고 같은 버전의 cupy-cuda12x를 사용함. 호환 STAR2.7.2a와 BWA0.7.17은 `tools/gpu/`에 분리 설치하고 source SHA256·실제 binary SHA256을 install_manifest.json에 기록함.

GPU 사용 근거는 확인한 범위대로 구분함.

| 근거 | 확인한 실제 사용 | 이 전처리에 적용할 때의 범위 |
|---|---|---|
| [Mansueto et al., GigaByte, 2024: CannSeek](https://doi.org/10.46471/gigabyte.135) | Methods에서 Parabricks `rna_fq2bam`을 RNA-seq benchmark에 실행했다고 명시함. 본 데이터의 RNA variant calling은 CPU STAR/GATK 경로를 사용함. | GPU RNA 정렬의 독립 실험 사용 근거임. 인간 개인 HAP·현재 GTEx 옵션과의 동일성 근거는 아님. |
| [Samarakoon et al., Bioinformatics Advances, 2025](https://doi.org/10.1093/bioadv/vbaf085) | WGS10개로 Parabricks4.3.0-1 FQ2BAM mapping·sorting·duplicate marking·BQSR 및 variant calling을 CPU/DRAGEN과 비교함. L4·A100·H100 GPU 사용함. | 독립적인 GPU WGS 정렬/성능 평가 근거임. 현재 개인 HAP·HiFi·B200의 가속배수를 뜻하지 않음. |
| [Spirito et al., npj Genomic Medicine, 2026](https://www.nature.com/articles/s41525-025-00547-8) | 100가족298명의 short-read WGS에 Parabricks3.8 BWA-MEM으로 hg38 mapping 후 HaplotypeCaller 사용함. Methods의 Identification and prioritization of short variants에 명시됨. | 독립 연구의 실제 분석 활용 근거임. GPU STAR·커스텀 coverage/Zarr 사용 근거로 확대하지 않음. |
| [Zhu et al., Parabricks 개발 프리프린트, 2025](https://doi.org/10.1101/2025.07.23.666378) | GPU 기반 short/long-read·transcriptome 소프트웨어 suite를 소개함. | 개발진의 도구/성능 보고로 구분함. 독립 사용 논문과 별개임. |

- 공식 도구 계약: [GPU minimap2](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/minimap2), [FQ2BAM](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/fq2bam), [GPU STAR](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/rna_fq2bam), [호환 버전](https://docs.nvidia.com/clara/parabricks/about-parabricks/software-overview/output-accuracy-and-compatible-cpu-software-versions), [CuPy 설치](https://docs.cupy.dev/en/stable/install.html) 기준으로 구현함. GPU minimap2의 preset·anchor tie 처리 등 때문에 CPU 결과와 bit 단위 일치를 가정하지 않음.


### 준비 BAM 성공 캐시

- Illumina paired 입력의 route→name collate→입력 metrics 검증이 성공하면 BAM 옆에 `.prepared.json` 증명을 저장함. 기본 all_to_both 및 해당 paired 입력 경로에 적용함. 별도 strict hap_specific 준비 경로는 기존 동작을 유지함.
- 원본·index·decode reference·도구·구현·염색체 순서·관련 정책과 출력 지문이 같고 BAM quickcheck가 통과할 때만 재사용함. 큰 파일의 지문은 inode·크기·mtime·ctime이며 모든 내용의 SHA256 재검사와 구분함.
- 증명이 없는 과거 BAM, 입력/정책 변경, 손상·부분 출력은 다시 준비함. 새 성공 증명이 만들어진 뒤의 실패 재개부터 재사용 가능함.
- `MFM_PREPARED_BAM_CACHE=0` 또는 명시적 overwrite는 재생성을 요청함. 재사용 여부를 `[PREPARED-CACHE]` 로그와 `prepared_cache` 메타데이터에 남기며 과거 입력 집계 시간을 이번 실행시간으로 표시하지 않음.
- 두 HAP 성공 후에는 기존 수명대로 준비 BAM과 증명을 함께 정리함. 최종 모델 입력/QC/로그와 소유 원본의 정리 조건은 바꾸지 않음.
- 새로운 설치 의존성은 없음. 기존 `requirements.txt`와 conda yml을 그대로 사용함.

이하 날짜가 붙은 운영 기록은 **당시 실험·설정 기준**임. 예전의 고정4명·500GB 디스크 여유·대형 큐 검증 대기 설명은 현재 공용 local 큐의 기본값이 아님. 현재 실행 동작은 이 절을 기준으로 확인함. 서버 전용 MedAI 운영 큐가 이미 실행 중이면 현재 설정과 상태 경로를 유지하며, 위 명령으로 같은 데이터를 중복 실행하지 않음.

3,212개 sample( WGS-only 2,763 / WGS–RNA pair 449 )을 개인 HAP 좌표의 WGS·RNA 모델 입력으로 변환한다. v2는 GENCODE v47 및 GTEx/TOPMed v11 RNA reference를 사용한다.

## 가장 간단한 실행

```bash
conda run -n jj_mfm bash -lc '
  cd /path/to/process_code/v2
  PIPELINE_MANIFEST=/path/to/preprocess_manifest.tsv \
  INPUT_ROOT=/path/to/input OUTPUT_ROOT=/path/to/output \
  ./run_pipeline.sh sample-list --manifest /path/to/preprocess_manifest.tsv \
    --sample-list /path/to/sample_list.txt --input-root /path/to/input \
    --output-root /path/to/output --sample-count all
'
```

- 한 sample: `./run_pipeline.sh sample <S> --manifest ... --input-root ... --output-root ...`
- WGS만: `./run_pipeline.sh wgs <S>`
- RNA만 재개: `./run_pipeline.sh rna <S>`
- 실패 sample은 `runner_local.sh`가 status/retry 목록에 기록하고 다른 sample은 계속 진행한다.
- `CLEAN_INTERMEDIATES_ON_SUCCESS=1`이어도 WGS Stage 23 또는 RNA Stage 25 final package 검증 PASS 전에는 cleanup하지 않는다. cleanup backend는 삭제 직전 PASS 보고서를 다시 확인하고, RNA는 현재 HAP variant TSV의 SHA-256을 QC·count meta·ASE QC와 재대조하며 stale 보고서를 거부한다.

## 비율 QC 정책

기본 `QC_RATIO_POLICY=record_only`는 WGS mapped/paired fraction과 RNA paired fraction을 QC JSON에만 기록하고, 수치가 0.90 미만이라는 이유로 sample을 중단하지 않는다. 0.90은 문헌의 보편적 합격선이 아닌 과거 운영 참고선이다. 과거 hard-fail 재현이 필요한 경우만 `QC_RATIO_POLICY=strict`를 명시한다. 파일 손상·decode 실패·sample/reference 불일치·R1/R2 불일치·0 pair·routing count 손실·최종 산출물 무결성 오류는 계속 해당 sample을 차단한다.

## 핵심 계약

- Stage 01은 BAM/uBAM/CRAM/FASTQ를 canonical BAM으로 표준화한다. Illumina는 BWA-MEM2와 Picard MarkDuplicates(mark only), long-read는 minimap2를 사용한다.
- Stage 02는 reference·VCF·CRAM/BAM 22 autosome 계약과 RNA input(quickcheck, SM, pair fraction)을 함께 검사한다.
- WGS는 `09 개인 FASTA → 10 HAP 재정렬 → 11 1 Mb windows → 12 Compact Zarr` 순서다.
- RNA Stage 18은 한 번의 depth 계산으로 HAP별 전장 `coverage.zarr`와 WGS `kept=1` 구간에 1:1 대응하는 `coverage_windows.zarr`를 함께 남긴다.
- RNA Stage 20은 HAP FASTA 생성 VCF의 exon ±1 kb variant를 WGS read로 검증하고, Stage 21은 validated variant를 덮는 fragment만 ASE 근거로 센다. Stage 23은 schema v5 provenance를 요구한다.
- sample별 command·상태·시간은 `metadata/execution_provenance.jsonl`, `metadata/stage_runtime.tsv`, `_run_status/*.stages.tsv`에 기록한다.

자세한 단계별 입출력: [PIPELINE.md](../PIPELINE.md)  
외부 서버 입력·출력·shard 전달: [docs/실행과_출력구조.md](../docs/실행과_출력구조.md)
QC 기준의 근거·보수성·RNA pair fraction 정의: [docs/QC_근거와_보수성.md](../docs/QC_근거와_보수성.md)

## RNA 성능 옵션 (2026-10-01)

Stage18은 contig별 병렬 전장 depth 계산,
Stage20은 contig별 WGS variant 검증, Stage22는 HAP×contig별 병렬 bedcov를 사용한다. Raw 배열·CNV 통계·필터는 기존과 같다.
`RNA_COVERAGE_JOBS=auto`, `RNA_CNV_JOBS=auto`가 기본이며 샘플 CPU 예산과 실제
affinity를 넘지 않는다. 14코어/샘플에서는 각각 최대6/13 worker다. `=1`로 직렬 비교할 수 있다.
Stage22의 중앙값은 여전히 HAP 전체 유전자에서 계산한다. STAR·MarkDuplicates는 바꾸지 않았고, ASE count만 validated-HAP evidence 규칙으로 강화했다.

**coverage 정책:** Stage 18은 `pipeline1_update` 및 GTEx `bam2coverage.py`와 맞추기 위해
`samtools depth -g 1024`로 duplicate를 포함한다. ASE count와 CNV depth는 QC 목적상 duplicate를
계속 제외한다. 구형 duplicate-excluded coverage는 fingerprint/metadata 검증에서 거부되어 Stage 18만
다시 계산되며 STAR index·정렬 캐시는 alignment-only fingerprint로 재사용한다.
세부 구현·동일성 검증·실측 시간·로그: [docs/pilot5_복구가이드.md](../docs/pilot5_복구가이드.md).

## WGS Stage12 독립 검증 가속 (2026-10-01)

trace의 후보 window 이진 탐색, bounded endpoint batch 기반 depth 누적합,
supplementary depth pass의 중복 segment 검증 제거, SQLite WITHOUT ROWID와
worker당64MiB cache를 적용했다. 2026-10-20 기본 경로는 아래의 분할 정렬 방식으로 갱신됐으며 SQLite는 fallback으로 보존함. 원본 CRAM→read/segment SHA-256 multiset과
primary/supplementary depth 대조·변조 탐지는 계속 수행한다. 필터·저장 배열·압축은 그대로다.
task/하위 phase 시간·완료 수가 sample 로그에 즉시 기록된다.

기존 `TRACE_JOBS=auto`와 sample CPU 예산을 사용할 수 있고,
`MFM_TRACE_CPU_RESERVE`·`MFM_TRACE_BLOSC_THREADS`로 조절한다. 명시한 worker를
몰래 늘리거나 실행 중 aligner/validator를 재시작하지 않는다.
`MFM_TRACE_SQLITE_CACHE_MIB=64`는 worker별 임시 page cache이며 DB 전체를 RAM에 넣는 설정이 아니다.
실측 범위와 남은 병목: [docs/Stage12_검증가속.md](../docs/Stage12_검증가속.md).


## 2026-10-06 병합

공통 전처리는 캐시 fingerprint, RNA 선행 quickcheck, 동적 CPU 배분, catalog 최종 PASS 계약을 반영했다. 원본의 phase/HAP 자산/coverage/ASE hash 검증을 두 코드에 공유한다. original의 de novo·assembly·Clair3·기존 실행 기능과 cleanup backend는 보존한다.

저수준 original의 일반 기본값은 `CLEAN_INTERMEDIATES_ON_SUCCESS=0`임. 현재 서버의 MedAI 기본 config는 모델 전용 정리를 명시하여 full QC·source trace·정리 후 검증에 통과한 샘플의 소유 원본/중간 산출물을 제거함. 공용 reference는 보존함. 배치 진입점은 `bash run_handoff.sh run --lab medai --batch 001 --config 설정.json`이다. 전달용은 `bash run_pipeline.sh`를 사용하며 기존 모델 입력/QC/metadata/log 보존 정리 정책을 유지한다. WGS 개인 FASTA/재정렬/window/compact 단계는 original09~12, 전달용07~10이고 RNA는 original13~24, 전달용11~22이다.


## 서버의 중앙 데이터 경로
- 이 서버는 `/NHNHOME/korea/mfm_data/`를 데이터 관리 루트로 사용한다. 전체 분류·실제 위치는 해당 루트의 `metadata/storage/데이터_경로목록.tsv`와 `README.md`에서 확인한다.
- `preprocessed_data/final/medai_001/`은 현재 MedAI 출력, `raw/medai/`는 원본·캐시, `reference/training/`은 학습 FASTA, `reference/preprocessing/`은 현재 HAP·RNA reference, `reference/intervals/borzoi/GRCh38/`은 Borzoi BED이다.
- 실행 중 MedAI 저장소와 동봉 자산은 실제 위치를 유지하고 중앙 경로에서 연결한다. 원본 소유 기록·검증 이력과 reference 서열은 변경하지 않았다. `reference/profiles/`는 패키지별 기존 reference 역할을 보존한다.
- 기본 `./data`, `./output`, `./assets`는 중앙 카탈로그가 있는 이 서버에서 중앙 경로를 사용한다. 사용자 지정 경로는 우선한다. 외부로 복사한 전달본은 상대 경로를 유지하며 `MFM_DATA_ROOT`로 별도 중앙 루트를 지정할 수 있다. 상태·잠금은 기존 run_state 위치를 유지한다.


## 내부 N·blacklist 보존 기본 정책 (2026-10-20)

- 정책 식별자: `annotated_acgtn_v1`. 염색체/HAP의 양끝 연속 N을 trim하고 원래 offset을 유지한다.
- `WINDOW_DROP_N=0`, `WINDOW_MAX_BLACKLIST_FRAC=1.0`: 내부 N과 blacklist가 있는 창을 보존하며 실제 구간에 ACGT가 없는 창만 저장에서 제외한다. window1MiB·seed42·기존 비중첩 좌표 정책은 유지한다.
- `metadata/wgs_windows.tsv`에 전체 후보와 제외 이유·N개수/비율/최대연속길이·blacklist비율·좌표대응률을 남긴다. 최종 compact의 windows.tsv에도 보존 창의 동일 주석이 포함된다.
- `metadata/reference_mapping/HAP*/chr*.chain.gz`와 `metadata/window_policy.json`은 최종 입력의 좌표 주석이다. 중간 FASTA를 영구 저장하지 않아도 개인→reference 위치를 재계산할 수 있다.
- window별 N 개수 `n_count`·비율 `n_fraction`·최대 연속 길이 `n_max_run`을 보존함. 비율의 분모는 padding P를 제외한 실제 길이임. N 존재/1% 초과/5% 이하 flag와 고정 학습 cutoff는 저장하지 않으며, 학습 소비자가 원값에서 선택 기준을 계산함. 내부 검사용 `n_nonacgt`도 유지함. 공통 입력 모듈에 N/P 정답 제외·유효 정답수 정규화를 구현함. 지정 모델의 학습 진입점 연결은 별도 확인 대상임. RNA/WGS 관측값은 N/blacklist 때문에0으로 바꾸지 않음.
- 기존 명시적 N0/blacklist 임계값 실행은 옵션으로 가능하며 구형 결과는 보존한다. 새 window schema와 입력 지문으로 호환성을 확인한다. source trace·최종 매핑/후보 원장 검증이 통과해야 원본·중간 파일을 정리한다.

## WGS 입력 통계 병렬 집계

- 대형 name-collated BAM은 samtools SAM 스트림을 QNAME 경계로 나누고 별도 프로세스에서 집계함. 첫 R1/R2·singleton·other·HP/XU/OD 판정은 기존과 동일함. 새 SAM·분할 BAM 파일은 저장하지 않음.
- 실제 CPU grant·affinity·Slurm 할당과 RAM 예약을 우선함. 72CPU·decode4 기준 생산자2+worker66까지 사용 가능함. 작은 입력은 작업 수와 실행 비용에 맞춰 직렬 처리함.
- 기본 chunk8MiB·단일 QNAME64MiB 상한과 worker당 대기2건으로 RAM을 제한함. 실행 전 RAM 피크를 계산하고10% 여유·worker 주소공간 제한·기존 RSS 보호를 적용함. 입력이 상한을 초과하거나 도구/worker가 실패하면 부분 집계를 완료값으로 저장하지 않음.
- routing QC의 input_metrics_execution과 INPUT-METRICS 로그에 예산·worker 요청/사용 수·wall/CPU시간·worker RSS를 기록함. 이는 집계 숫자와 구분된 실행 정보임.


## 전과정 검증·재개·속도 보완 (2026-10-20)

| 항목 | 현재 구현 |
|---|---|
| 소형 E2E | 실제 BAM/CRAM·FASTA를 compact→주석/매핑→독립 source trace→최종 QC→정리 후 재조회로 검증함. upstream 정렬·variant calling은 준비된 fixture를 사용하며 전부 실행한 검사와 구분함. |
| 안전 재개 | compact·최종 주석·contig/HAP별 source trace의 성공 증명을 각각 저장함. 입력·설정·구현·출력 지문이 일치한 성공 결과만 재사용함. 주석 실패 후 compact 재생성과 완료 HAP 재검증을 생략함. |
| 개인 VCF | `VCF_EXTRACTION_MODE=chromosome`을 기본 사용함. 전체 다인 panel concat 전에 염색체별 개인 non-reference를 임시 BCF로 추출하고 seq/tag를 생성함. 개인 결과만 합침. `concat`은 기존 방식 선택임. |
| 4bit packing | 4MiB/65,536 read 제한의 일괄 연산 사용함. 홀수 길이 padding·offset·지원 문자·오류 위치를 유지함. HG00189 실제 read20만 개/CPU1개에서 중앙값2.866→0.140초, 약20.4배임. 전체 전처리 가속률과 구분함. 실측 최대 RSS261MB임. |
| 학습 공통 입력 | `model_input.load_window`로 packed DNA·WGS·RNA 좌표와 실제 길이를 대조함. N/P DNA 정답 제외·유효 정답 수 정규화 함수 제공함. 특정 모델의 학습 진입점에 연결한 상태와 구분함. 원 WGS/RNA depth는 보존함. |
| 관측 깊이 QC | 신규 WGS는 `compact/<layout>/chr*/window_quality.tsv`에 depth 합계·평균·0 coverage 비율·최댓값 기록함. 신규 RNA는 `coverage_windows.zarr`의 `depth_mean/depth_zero_fraction/depth_max` 배열에 기록함. padding 제외한 실제 길이로 계산하며 저장 필터로 사용하지 않음. |

완료 증명은 `metadata/pipeline_steps/{compact,annotations}.json`, 독립 trace는 `qc/stage12_trace/<hap_read_mode>/hap_steps/chr*/HAP*.json`에 보존함. 큰 입력/출력은 inode·크기·수정/상태 시각 지문을 사용하며, 작은 입력은1MiB 이하 내용 SHA256도 대조함. 모든 Zarr chunk를 재해시하는 검사와 구분함. 최초 독립 source trace·최종 QC는 유지함. 출력 오염/검증 실패는 성공으로 채택하지 않음.

VCF 성공 캐시는 원 panel·index 지문, 정책, 검사 구현, 압축 VCF/index 전체 SHA256이 일치해야 사용함. 기존 무증명 VCF는 실제 검사를 다시 수행하되 생성 당시 원 panel까지 입증한 것으로 기록하지 않음. 입력/정책이 달라지면 명시적 재생성이 필요함. 임시 개인 BCF·염색체별 VCF는 호출 종료 시 정리함. 전체 CPU 풀은 유지하며 동시 개인/염색체 작업을 실제 CPU grant와 RAM 예산 안에서 나눔. 예약 추정과 cgroup 강제 제한은 구분함.

공통 로더는 최종 QC PASS와 full compact만 기본 허용함. RNA는 배열 순번 대신 contig·개인 시작/끝 좌표로 연결함. 입력 코드의 N=4/P=5는 DNA 정답에서 제외함. `torch_dna_cross_entropy`는 유효 정답 수로 정규화하며 전부 제외된 batch도 gradient가0인 loss를 반환함. PyTorch는 학습 환경의 선택 의존성이며 전처리 requirements에 추가 설치하지 않음. RNA/WGS loss 선택은 해당 모델 목적에 맞게 별도 결정함.

```bash
# 완료 패키지의 모델 공통 입력 검사
python tools/check_model_input.py /완료/샘플 --contig chr22 --hap HAP1 --index 0
# 과거 완료본은 수정하지 않고 window 관측 깊이만 별도 집계함
python tools/summarize_windows.py /완료/샘플/compact/LAYOUT --contigs chr22 --output /별도경로/깊이.tsv
# 실제 엔진·중단/재개·최종 QC·정리 회귀 검사
PYTHONPATH=src:vendor/MFM-DNARNAFM:. python -m unittest discover -s tests -p test_pipeline_e2e.py
```

구완료본에 새 sidecar/배열/완료 증명을 소급 추가하지 않음. 현재 실행 중 프로세스에 새 코드가 소급 적용된 것으로 간주하지 않음. HG00189는20:41 RNA·최종QC/정리·용량실측 완료함. 새 정책 HiFi 대형 실행, 지정 모델의 실제 training step, 외부 서버 신규 환경 설치·Slurm 실측은 후속 운영 확인 대상임.


## 후속 QC·반복 입력 조회 최적화 (2026-10-20 21시)

- **출력 BAM/CRAM QC를 양쪽 공통 병렬 모듈로 통합함.** `realign/output_scan.py`에서 index가 있으면 contig별로 나누고 집계함. 미배치 `*` 전체를 포함하며 여러 미배치 CRAM container는 임시 CRAI만 정규화함. 원본 CRAM/CRAI는 보존함. index가 없거나 자동 선택 시32MiB 미만이면 순차 전체 scan함.
- CPU grant·affinity·Slurm과 동시 sample/HAP 수를 반영함. RAM은 최대1TB·예약900GB와 더 작은 할당을 우선함. 입력 크기·최대 reference 길이·worker buffer·10% 여유로 피크를 계산함. worker 주소 공간 제한과 실제 RSS 관측을 사용하며 cgroup hard limit으로 표현하지 않음.
- `realign.json`의 `output_qc_execution`과 OUTPUT-QC 로그에 실제worker·예산·wall/CPU시간·worker RSS를 기록함. QC 숫자와 구분함. 다중 미배치 CRAM·BAM/CRAM 직렬/병렬·samtools flagstat 일치 검사 통과함. 전체 WGS 대형 CRAM의 실제 속도 배수는 다음 실행 측정 대상임.

반복 학습 조회에는 아래 reader를 사용함. 기존 `model_input.load_window`는 단건 조회 API로 유지함. 특정 모델의 학습 진입점은 아직 지정되지 않았으며 이 reader 제공과 별도로 연결해야 함.

```python
from integrate_wgs_rna.window_reader import WindowReader

reader = WindowReader('/완료/샘플', max_cached_contigs=2,
                      max_metadata_bytes=64 << 20,
                      max_data_cache_bytes=64 << 20)
x = reader.load('chr22', 'HAP1', 0)  # 저장 배열 순번임.
# x: dna_codes, dna_target_mask, observation_mask, wgs_depth,
#    actual_length, raw_start, rna_depth
print(reader.cache_info())
reader.close()
```

- 완료 후 불변인 local Zarr 패키지에 사용함. QC PASS/full manifest 확인 후 좌표·Zarr array metadata를 제한적으로 재사용함. 해제된 chunk는 별도 LRU의 최대64MiB까지만 재사용함. 전체 DNA/depth 배열을 상주시키지 않음. `max_data_cache_bytes=0`이면 chunk 재사용을 끔. 출력 배열은 독립된 값이며 캐시를 수정하지 않음.
- worker마다 reader를 사용함. pickle/spawn에는 경로·설정만 넘기며 fork 후에도 worker에서 새로 연결함. QC/manifest·좌표·캐시 chunk의stat변경은 실패함. stat지문 검사는 모든 파일의 내용SHA 확인과 구분함. 외부에서 완료본을 수정하면서 reader를 사용하는 방식은 지원하지 않음.
- RNA는 배열 순번이 아닌 contig·개인 시작/끝 좌표로 연결함. N/P DNA 정답 mask·실제 길이·원 WGS/RNA depth를 유지함. LRU상한·중복 RNA 좌표·길이오염·fork/spawn 재연결·단건API대조 검사 통과함.
- 실제 HG00189 chr22/HAP1의 연속16창을3회 비교함. reader 연결 포함 중앙값0.570→0.127초(약4.48배), 모든 반환 배열 일치함. CPU1·RSS381MB임. 연속 조회의 chunk재사용 효과이며 무작위 전장 접근/전체학습/전체전처리 속도 배수로 사용하지 않음.

완료 샘플을 전수 조회하려면 다음 명령을 사용함. 검사 결과는 완료 패키지 바깥에 저장함.

```bash
python tools/audit_model_inputs.py /완료/샘플 --output /별도경로/입력검증.json
# --workers 생략 시 실제 CPU/RAM/44개 HAP축에서 최대 병렬도 선택함.
# --max-per-axis 3은 축별 표본 검사이며 기본값0은 전수 검사임.
```

HG00189는20:41에 독립trace44개·RNA·full package105검사·정리후검증을 통과함. 별도공통reader에서도 저장5271창 전수검사 PASS임. 실제5,505,152,296bp·N37,076,764bp·padding21,891,800bp를 읽어 원장과 대조함. N포함360창·padding포함44창·RNA5271창의 대응을 확인함. 이 중 N은 보존된 서열이며 DNA 정답에서만 제외함. 전수검사는 reference/BAM을 다시 읽는 source trace와 구분함.

실물 용량은 동일한 `compact+rna_personal` 일반파일 합계 기준29.263316403GB→32.365483201GB, +3.102166798GB(+10.600872%)임. RNA내부QC/JSON을 포함하고 sample 최상위metadata/qc/logs는 제외함. cleanup의 엄격한 모델선택 부분집합은32.364783409GB임. 한 명의 관측치이며 동결118.215TB·배치 배정이나 전체 모집단 추정에 일괄 곱하지 않음. 새 정책 HiFi 전장·대형 CRAM 실측·지정모델 training step·외부 신규환경/Slurm 검증은 후속 항목임.


## RNA HAP 병렬 처리·독립 source trace 가속 (2026-10-20)

RNA Picard 중복 표식을 HAP1·HAP2에 동시에 요청함. 기본 `RNA_MARKDUP_HAP_JOBS=2`이며 `1`이면 직렬 실행함. 두 HAP가 모두 성공한 뒤 다음 단계와 `RNA_STOP_AFTER`를 처리함. 실패한 경우 이미 성공한 HAP의 유효 캐시는 보존하며 다음 단계는 차단함. 부모만 sample/stage 완료 상태를 갱신하고 HAP별 명령·로그·실행 시간을 분리 기록함. 종료 신호를 받으면 해당 자식과 CPU/RAM lease를 정리함.

Picard 한 작업의 기본 CPU 요청은 최대4개임. Java와 후속 BAM 인덱스·집계가 실제 grant를 사용하도록 수정함. 두 HAP를 실행한다고 전체72 CPU 풀을 축소하지 않으며 남는 자원은 공유 풀의 다른 작업이 사용함. Picard 자체의 직렬 계산 구간 때문에 CPU grant 전부가 계속 계산에 사용된다는 의미는 아님. heap 기본값은 기존 메모리 profile에 따르며, 작업별 예약은 `(heap GiB + native/buffer 2GiB) × 1.1`로 계산함. heap8GiB이면 약11.81GB/작업, 두 작업 합계 약23.62GB임. 전체 RAM 상한1TB·예약900GB와 더 작은 Slurm/lease 한도를 우선함. 실행 계획·입력 BAM 크기·wall time을 `hap*.markdup.json`의 `resource_execution`에 기록함.

독립 trace는 저장된 read의 4bit 해제를 배치 처리하고 기존 canonical SHA256 바이트를 그대로 계산함. source 측 투영·서열 계산과 분리된 구현을 유지함. 기존 SQLite 랜덤 갱신 대신 SHA256 원값을 순차 기록한 뒤 한도 내에서 분할 정렬해 정체성과 중복 개수를 모두 대조함. 누락·추가·서열/metadata/CIGAR/segment/depth 오염 검사와 전수 대조를 유지함.

| 설정 | 기본값·동작 |
|---|---|
| `MFM_TRACE_MULTISET_BACKEND` | `auto`: 예상 임시 공간이 가용량의90% 안이면 `partitioned`, 부족하면 기존 `sqlite`를 선택함. 둘 다 부족하면 실행 전에 중단함. 명시값 `partitioned`·`sqlite`도 지원함. |
| `MFM_TRACE_SORT_LEAF_MIB` | 64MiB 원자료 leaf 한도임. 불일치 집계 배열을 포함해 최대8배 공간을 작업 RAM 계획에 추가함. 큰 leaf는 다음 SHA byte로 재분할함. |
| `MFM_TRACE_DECODED_CACHE_MIB` | 작업당 Zarr 해제 chunk 캐시64MiB임. |
| `MFM_TRACE_PREFETCH_PACKED_MIB` | packed read 배치 최대16MiB임. 긴 read에서는 배치 read 수를 줄임. 단일 read 자체가 한도보다 크면 명시적으로 중단함. |
| `MFM_TRACE_SQLITE_CACHE_MIB` | SQLite fallback의 기존 작업별64MiB cache임. |

`memory_preflight_trace.json`에 read·segment·CIGAR 수 기반 scratch 예상량, filesystem 종류, 선택 backend와 worker RAM 한도를 기록함. tmpfs scratch는 RAM 예약에 임시 파일 예상량도 포함함. 다른 작업의 디스크 공간까지 전역 예약하는 기능은 아니므로 동시 디스크 사용량에 따라 실제 여유가 달라질 수 있음. 자원 예약·주소 공간 제한·RSS 보호는 cgroup 전체 RAM 강제 제한과 구분함.

실제 HG00189 chr22/HAP1 read20만 개에서 기존·수정 SHA256 전체가 일치함. read 해제·digest와 multiset 비교 시간을 각각 측정했으며, 전체 염색체 source trace 시간과 구분함. 실제 소형 Picard BAM에서도 직렬/동시 실행의 모든 alignment와 duplicate 표식·metrics가 일치함. 성공·부분 실패·START/STOP·종료 신호·lease 반환과 분할 정렬의 반복 digest/불일치 예시/긴 read 메모리 경계 검사를 통과함. 대량 실행 재개와 완료 데이터 재작성은 수행하지 않음.

실물 표본을 CPU1개에서 각3회 측정한 중앙값은 read 해제·digest 1.448→0.635초(2.28배), multiset 비교 1.201→0.320초(3.75배)임. leaf1MiB로 강제 분할해도 0.446초(2.69배)이며64개 leaf로 실제 분할됨. 프로세스 누적 peak RSS는 171.3MB임. 전장 source trace의 wall time 실측과 구분함.


## 2026-10-20 후속 실행·STAR/ASE 최적화

ASE 입력에서 반복 reference 이름 생성을 줄이고, 순환 참조가 없는 read table 구축 중 cyclic GC만 일시 중지함. 호출 전 GC 상태는 성공/예외 모두 복원함. 기존 AS·MAPQ·duplicate·mate 선택 순서를 유지함. STAR 길이 보정은 같은 수식을 두 mate 직접 비교로 계산함. 실물 참조 RNA prefix100만 alignment의 두 입력 읽기에서 기존·수정 전체 table의 canonical digest가 같았고, 각3회 중앙값8.750→4.148초(2.11배)임. chr1 첫50만 alignment로 만든208,255 fragment의 길이 보정은0.246→0.129초(1.90배)이며 모든 보정값이 같음. 전체 개인 HAP ASE 시간 가속률과 구분함.

| 설정 | 동작 |
|---|---|
| `RNA_STAR_INDEX_HAP_JOBS` | 일반 기본2임. 1이면 번호별 경로로 직렬 실행함. 두 HAP index를 공유 CPU/RAM 풀에 동시 요청함. 실제 grant와 단계별 예약을 적용하며 둘 모두 완료 후 STOP/다음 단계를 처리함. |
| `RNA_STAR_INDEX_BENCHMARK` | 기본0임. 1이면 같은 개인 HAP/GTF로 직렬·동시 index를 모두 생성하고 주요 파일 SHA256을 비교함. 동일하고5% 이상 빠르면2를 추천함. 그 외 기존 직렬 캐시를 복원함. 결과는 `qc/rna_index_benchmark.json`에 보존함. 추가 생성/해시 시간은 일반 처리시간과 구분함. |
| `RNA_ALLELIC_JOBS` | 일반 기본auto임. 양수는 상한으로 사용함. 큰 입력은 아래2026-10-31의 QNAME shard 경로를 사용하며 legacy에서 Linux fork table 공유를 사용함. CPU/RAM 범위에서 실제 동시 수를 정함. |
| `RNA_ALLELIC_BATCH_SIZE` | 기본50,000 fragment/task임. 작은 입력·단일 CPU·안전한 fork 조건이 아니면 직렬 처리함. |
| `RNA_ALLELIC_BENCHMARK` | 기본0임. 1이면 입력 table과 fragment count/stat 전체를 순차 계산 결과와 대조하며 소요 시간을 기록함. 전체 비교 비용은 일반 처리시간과 구분함. |

fragment 하나의 모든 mate와 두 HAP 증거는 같은 key 작업에서 계산함. 두 HAP duplicate 이름의 전체 합집합을 먼저 제거하고 기존 AS 비교·validated variant·공통 gene·split-half 규칙을 사용함. count는 정수로 합산하고 beta-binomial/BH FDR는 기존 최종 단계에서 수행함. 대형 read dict의 pickle 전달은 사용하지 않음. 2026-10-31부터 새 네이티브 binary QNAME shard를 짧은 임시공간으로 사용하고 자동 제거함. 기존 Python table 직렬 경로도 명시 선택 가능함.

fork worker의 주소 공간과 공유 페이지를 포함한 RSS를 보수적으로 예약함. 유휴 worker도 종료까지 RAM 예약을 유지하며, 실제 CPU 사용량과 CPU grant를 구분함. 기존 pool audit에 `actual_cpu_seconds`, `average_used_cores`, `sampled_leased_cpu_seconds`, `cpu_fraction_of_lease`, `peak_tree_rss_observed_bytes`, `peak_single_process_rss_bytes`를 추가함. CPU는 wait한 자식/후손의 누계이며 tree RSS는 RAM guard의 실제 시각 기반 표본이며 잠금 대기로 관측 간격이 늘어날 수 있음. 순간 피크·cgroup 전체 강제 상한과 구분함.

HG00275를 전체 실물 검증 대상으로 선정함. 원본 다운로드가 끝나면 최신 코드로 WGS→RNA→full/source trace→정리 후 검증을 실행함. 같은 샘플에서 STAR 직렬/병렬과 ASE 동등성을 비교한 뒤, 유효하고5% 이상 빠른 병렬 방식만 후속 MedAI001 신규 큐에서 선택함. 보존 archive28명은 이번 신규 큐에서 제외함. 전체 CPU72·RAM1TB/예약900GB를 유지하며 더 작은 Slurm 할당을 우선함. 후속 큐는 입력/중간 파일 예상 공간과500GB 여유로 입장을 제한함. 현재 전장 실측과 후속 대량 큐는 실물 검증 완료 대기임.


### 2026-10-21 입력 준비·사전 검증 최적화
- Stage01 전수 flag·HP/PS 통계를 인덱스 없이 samtools flagstat와 tag 존재 count로 계산함. QC fail·unmapped·secondary/supplementary·duplicate·HP=0/PS=0도 포함하며 원본/최종 BAM 전체 카운트를 대조함. 64MiB 미만은 기존 Python 경로를 사용함. `MFM_RECORD_STATS_BACKEND=auto|python|native`로 선택 가능함.
- HG00189 WGS 앞100만 레코드의 같은 BAM 표본에서 각3회 중앙값1.435→0.099초(14.5배), 모든 카운트 동일함. 이 수치는 집계 부분이며 전체 Stage01 가속률이 아님.
- CRAM view 디코딩과 sort의 추가 thread·각 main thread 합을 실제 grant 이내로 배분함. 전수 집계는 기존 Stage01 임대 범위와 RAM 예약을 유지해 후속 sort/index가 자원 없이 실행되지 않게 함.
- Stage02 reference M5는 case를 보존한8MiB chunk 해시와 chromosome 병렬 계산을 사용함. VCF REF 검사는 genotype를 제거한 모든 site/allele에 기존 `norm --check-ref e`를 수행함. 해당 샘플 GT·위상 검증은 Stage03 이후 기존 계약을 유지함.
- `PREFLIGHT_JOBS`는 CPU 풀 안에서 검사 동시성을 지정함. VCF pipeline당 최대2CPU·RAM2GiB, M5당256MiB와 부모4GiB를 보수적으로 예약함. 공통 성공 cache 기본 경로는 `$DATASET_ROOT/_validation_cache/input_preflight`이며 `PREFLIGHT_CACHE_ROOT`로 변경 가능함. 파일·index·reference·도구의 inode/크기/nanosecond mtime·ctime와 검사 옵션이 같을 때만 재사용함. lock·원자적 저장·검사 중 입력 변경 검출을 적용함.
- Stage02도 throughput 모드에서 전체72CPU 풀을 사용할 수 있음. 실제 활성 작업 수·완료에 따라 CPU를 반환하며 생성된 작업의 RAM은 종료까지 예약함. RAM1TB/예약900GB·더 작은 Slurm 할당 우선 원칙을 유지함.


### STAR 두 HAP의 실제 동시 실행 확인 (2026-10-21)
STAR index 직렬/병렬 비교는 같은 총 CPU 예산을 사용함.72CPU가 가용하면 직렬 HAP당 최대72, 병렬 HAP당 최대36으로 요청하여 첫 작업의 전체 풀 선점을 방지함. 더 작은 명시 예산도 두 작업에 나눠 적용하며 함수 호출 후 기존 상한은 원복함. 전체 CPU풀72와RAM1TB/예약900GB는 유지함.

실제 STAR genomeGenerate 명령의 시작·종료 시각과 thread를 `hap*.index.json`의 `index_execution`에 기록함. 보고서의 `actual_overlap_seconds>0`, `cpu_budget_comparable=true`, 핵심 index SHA256 일치와 wall time5% 이상 개선을 모두 만족할 때만 병렬2를 추천함. 겹침 증거가 없거나 더 많은 CPU를 사용한 비교는 직렬1을 유지함. 일반 기본은1이며 HG00275 전장 실측에서 최종 선택함.


ASE fragment 병렬 계산은 worker의 PID·시작/종료·CPU시간으로 실제 동시성을 확인함. metadata에 최대 활성 worker·부모/자식 CPU 누계·평균 사용 코어를 기록함. 실제 두 worker 이상의 overlap이 없으면 후속 큐 선택용 `assignment_speedup`은1로 두며 raw 비교값은 assignment 하위에 보존함. 입력 읽기는 기존 직렬 최적화를 유지하고, Counter·통계·최종 beta-binomial/BH 계산은 바꾸지 않음.

## 재정렬 하위 단계 자원 배분 (2026-10-21)

- 기본 `MFM_REALIGN_PHASE_RESOURCES=auto`임. 단일 샘플 Illumina BWA·name-collated 입력에 적용함. `0`은 기존 방식, `1`도 적용 가능한 입력 조건에서만 활성화함. 다중 샘플 `--jobs>1`·long-read·dry-run은 기존 실행을 사용함.
- 입력 준비 뒤 같은 CPU lease를 index/정렬/마무리 구간에 맞춰 원자 갱신함. index는 최대2CPU를 예약하고 나머지는 다른 샘플에 반환함. 부모·자식 lease를 중첩하지 않음. 실행 중 BWA thread 수는 변경하지 않음.
- `MFM_BWA_INDEX_HAP_JOBS=auto|1|2`이며 기본 auto임. HAP별 reference 길이로 `40×bp+2GiB`에10% 여유를 둔 작업 RAM과 부모4GiB를 먼저 예약함. CPU/RAM이 허용하는1개 또는2개만 실행함. cache별 flock·완료 manifest 검증은 유지함. worker 주소 공간 제한과 전체 RSS 보호를 함께 사용함. cgroup 강제 상한과는 다름.
- 정렬 시작 직전 현재 가용 CPU/RAM에서 HAP 동시 수·BWA·sort thread를 원자적으로 재산정함. 기존 `MFM_WGS_*_REQUESTED=auto` 요청만 재계산하며 숫자로 명시한 요청과 `WGS_AUTO_THREADS=0`은 보존함. auto HAP은 CPU16 이상에서2개를 우선하며 부족하면1개로 시작함. 실행 중 BWA의 thread는 변경하지 않음.
- name-collated FASTQ 변환·duplicate 복원·sort 주 프로세스에 HAP당3CPU를 계상하고 sort 보조 thread까지 예약함. RAM은 HAP 동시 수·reference 길이·sort 버퍼로 계산하며 최소192GiB의 기존 보수적 예약을 유지함. 예산 부족 시1CPU/4GiB만 보유하고 대기함.
- 종료된 manager/child PID의 예약은 프로세스 그룹까지 확인한 뒤 회수함. 분리 세션에 실제 계산이 살아 있으면 예약을 보존함. 선택한 HAP 수·thread·예약 CPU/RAM·대기시간을 resource_phases/effective_threads 및 lease에 기록함.
- CPU 풀은 허용된 전체 범위를 사용함. RAM 상한1TB·예약900GB이며 Slurm 할당이 더 작으면 해당 CPU/RAM을 우선함. index1작업도 예산을 넘으면 실행 전에 실패함.
- `qc/<sample>.stage11_timing.json`에 `resource_phases`, `bwa_index_jobs`, `bwa_index_execution`, `effective_threads`를 기록함. `logs/personal_realign/*.bwa_index.resource.json`에 실제 index CPU·시간·RSS·주소 공간 제한을 기록함. cache 재사용은 새 생성 실측과 구분함.
- 작은 실제 BWA 입력의 이전/신규 read별 결과와 QC 동일성을 검사함. 전장 데이터에서의 처리시간 개선율은 별도 실측이 필요함. 추가 환경 패키지는 없음.


## 다운로드 전송률 로그

S3를 포함한 HTTP/FTP 다운로드는 실제 본문 수신 바이트와 요청 소요시간을 기록한다.
로그는 실행 상태의 `stats/download_transfers.jsonl`과 입력 저장소의
`.download_records/transfer_rates.jsonl`에 남는다. 시작·첫 구간·이후 약30초 간격·완료/실패를 기록하고
기존 명령 로그에도 `[download-rate]` 행으로 MB/s와 Gbps를 출력한다.

- `received_body_bytes`, `network_seconds`, `network_MBps`, `network_Gbps`: 연결·요청·본문 수신 기준. 헤더/TLS 바이트와 디스크 쓰기·재시도 대기·최종 검증은 제외한다. 재수신한 본문과 원격 식별/재개 확인용 본문도 실제 수신량에 포함한다.
- `resume_offset_bytes`, `new_file_bytes`: 기존 부분 파일과 이번에 새로 저장한 바이트를 분리한다.
- `transfer_seconds`, `useful_transfer_MBps`: 재개 검사·재시도 대기·로컬 쓰기까지 포함한 전송 단계의 실효 처리율이다.
- `verification_seconds`, `retry_count`, `status`, `failure_type`: 최종 크기/checksum·alignment 검사 시간과 재시도·실패 상태다.
- 검증된 기존 파일을 재사용하면 새 다운로드 속도를 생성하지 않는다. query·서명·사용자 정보는 새 전송률 로그에서 제거한다.

성공 영수증 `.download_records/<파일해시>.json`의 `transfer_metrics`에도 최종 통계를 보존한다.
새로 시작하는 다운로드 프로세스부터 적용되며 이미 실행 중인 프로세스의 코드는 교체하지 않는다.
동시 다운로드 전체 처리율은 개별 전송시간을 합하지 않고 같은 관측 구간의 총수신량/경과시간으로 계산해야 한다.

### 전체 샘플 GPU 대조의 검증 범위 (2026-10-29)

- HG00275(short-read WGS+RNA pair)의 전체 원본18.05GB와 chr1–22/HAP1·2를 사용함. 기존 CPU 완료본 보존, 별도 GPU 결과/상태/입력 경로 사용함. 다른 샘플·GPU 작업은 중단하지 않음.
- 전체 전처리·full QC/source trace 이후 HAP 서열/좌표/window 주석, WGS/RNA1bp depth, readset 수량과 단계별 실측을 대조함. readset은 정렬 차이/저장 순서 때문에 단순 바이트 동일성을 강제하지 않음.
- RNA는 같은 GPU alignment를 CPU samtools depth로 전수 재계산하여 실제 저장된 GPU depth의 SHA256과 대조함. 같은 STAR2.7.2a·Normal·SJ 삽입 한도·index overhang·FASTQ·thread 조건의 CPU RNA 정렬도 추가하여 도구/정책 차이를 분리함.
- 다운로드·실패 재시도·대조 검증 시간은 모델 전처리 시간과 구분함. 과거 CPU와 신규 GPU의 CPU grant/cache/실행 부하가 달라 기존 기록 대비 배율은 관측값으로 보고함. HAP 병렬 단계는 시간을 합산하지 않고 합집합으로 계산함.
- 이번 대조만 `http_range_workers=16`을 사용함. 기본1이며 허용1~16임. 각 Range의 버전/길이/XML 오류를 검증한 뒤 연속 prefix 순서로 기록하고 전수 CRC/공식 WGS MD5를 확인함. 최대16×32MiB 결과 및 수신 버퍼는 기존 다운로드4GB 예약 안에서 제한함.
- 전체 대조가 통과할 때까지 비교용 원본/중간 파일을 보존하고, 통과 후 기존 정리·정리 후 무결성 검사 및 소유권 확인 원본 삭제를 수행함. 전장 결과와 가속률은 실행 완료 후 판단함.

### 2026-10-29 CPU·RAM 운영 정책

- RAM 기본 상한은 1.5TB이며 10% 여유를 둔 예약 예산은 1.35TB입니다. 물리·Slurm·cgroup 할당이 작으면 우선합니다.
- CPU는 실제 할당 한도까지 사용합니다. RAM과 최소 CPU를 함께 확보할 수 있는 대기 작업에 분배하고, 자동 WGS 정렬은 최소 4개 계산 스레드를 확보합니다.
- 내장 자원 감시기가 공유 풀의 최신 승인 정책, 대기 작업 할당, RSS 감시를 유지합니다. 기존 계산 작업은 보존하며 승인된 HAP 경계에서 자원을 재배정합니다.
- 새 cgroup 권한이 없으면 예약·작업별 주소 공간 제한·RSS 감시로 보호하며, 이 방식은 전체 프로세스의 커널 강제 상한과 다릅니다.

### 실제 유휴 CPU를 활용하는 추가 투입

- maximum 모드에서 CPU가 예약만 되고 실제로 비어 있으면, 10초 관측과 최근 2초 평균 최대 사용량을 확인하여 고정 자원 단계에 추가 대여합니다.
- 실제 CPU cgroup quota와 풀 상한이 같을 때만 동작합니다. 대여 단계는 nice 19로 실행해 기본 작업에 우선권을 줍니다. RAM 예약·RSS·cgroup 여유는 추가 작업까지 합산합니다.
- 기본 CPU 예약은 풀 한도, 유휴 대여는 별도로 최대 한 풀입니다. 논리 스레드 합이 68을 넘을 수 있으나 현재 커널 CPU 상한은 68코어입니다. 예약과 실제 사용을 구분하세요.
- watchdog.json에 기본 예약·대여·실측 CPU 및 작업별 대기 사유를, cpu_backfill_events.jsonl에 추가 투입 근거를 남깁니다.

### 앞단 IO 선택과 전체 BAM 비교 (2026-10-29)

- 선택 인자: `--io-backend cpu|nvcomp` 사용 가능함. 기본 CPU이며 GPU 압축·CRAM CRC 선택은 별도 실험 경로임. GPU 정렬·coverage·QC·서열 packing 선택과 독립임.
- 공용 장치: 현재 비교 실행은 GPU `0,1,2,3` 및 동일 `gpu_pool_dir`를 사용함. 샘플별 장치 고정 없이 준비된 HAP가 빈 GPU를 예약함. CPU 선행 단계에서는 GPU가 유휴일 수 있음.
- GPU BAM 준비: CPU samtools comparator와 read·태그·순서를 보존하고 BGZF 압축만 nvCOMP6.0.0.12로 수행함. 한 writer의 CUDA context와 codec을 전체 스트림에서 재사용함. bounded 배치·입력 prefetch·GPU 압축·독립 CPU 해독/CRC 전수 대조를 겹쳐 수행함. GPU 요청은 최소6CPU를 확보하며 CPU로 조용히 대체하지 않음.
- 안전한 입력 검사: 원본 BAM의 압축해독/CRC와 MD5·VCF REF 검사는 CPU 경로를 유지함. CRAM CRC만 선택 GPU 배치 경로를 지원함. 손상된 원본을 nvCOMP로 직접 해독하여 검증하지 않음.
- 자원 대기: RAM·최소 CPU가 실제로 들어갈 수 있는 대기 요청의 우선권을 보존함. 들어갈 수 없는 요청 때문에 검증 작업의 CPU 재확보가 영구 차단되지 않도록 함. 고정 base RAM도 공유 예약 예산 안에서 확보함.
- 소규모 실측: 실제 HiFi 약1.9GB 준비된 BAM 스트림의 압축/쓰기 CPU4코어7.220초·GPU8.058초, HG00275 전체 원본CRAM14.807GB CRC CPU13.471초·GPU15.580초임. 전체 BAM 준비·전체 전처리 가속률로 확대 해석하지 않음.
- 전장 비교: NA18959 원본83,941,493,808바이트 전체를 동일 Stage01 함수로 GPU IO 준비하고 기존 CPU canonical BAM과 모든 record 바이트·태그·순서·reference header·CRC를 독립 대조함. CPU/GPU의 압축파일 바이트와 실행 @PG header는 달라도 됨. 동일6CPU 추가 CPU 대조도 수행함. 전처리 운영의 샘플 CPU 상한을6으로 바꾸는 설정이 아님.
- 시간 해석: 기존 NA18959 Stage01은 실제 계산967.711초(16분7.7초), CPU 예약 대기754.202초(12분34.2초)임. 할당19CPU·평균 실제7.2363코어임. 전장 GPU 실측 완료 전 가속률을 확정하지 않음.

### 2026-10-29 최신 RAM 완화·계산 우선 정책

- RAM 운영 상한을 1.65TB로, 10% 여유 후 예약 예산을 1.485TB로 변경했습니다. 물리·Slurm·cgroup 중 더 작은 실제 할당을 우선합니다. 현재 cgroup 상한은 1,717,986,918,400바이트입니다.
- STAR index 생성과 정렬 RAM을 분리합니다. 실제 index 파일·FASTA 크기, STAR 생성/정렬 상한, 스레드 버퍼에 15% 여유를 합산하고 최소 128GiB를 확보합니다. 기존에 일괄 예약하던 384GiB를 입력별로 정정합니다. 실행 중 RAM 증액은 거부하며 계산 인자와 완료 캐시는 유지합니다.
- 고정 heap MarkDuplicates도 유휴 CPU 대여를 사용할 수 있습니다. 기본 CPU 예약과 실제 사용량을 구분하고, 전처리와 같은 cgroup을 쓰는 별도 작업의 CPU도 따로 측정합니다.
- 현재 40TB 저장소 포화 재발을 막기 위해 공유 정책에 WGS 준비 단계 1개, 선행 다운로드 1개 및 다운로드 전 500GB 여유 조건을 적용했습니다. 이는 디스크 임시 산출물에 대한 운영 조건이며 CPU 풀 68개는 유지합니다.
- 제어기 중지 요청 뒤에는 신규 worker와 다운로드를 시작하지 않습니다. 기존 계산과 완료 원본·모델 입력은 보존합니다.

- BAM 준비 직전에는 실제 입력 크기×2+50GB와 최소350GB 디스크 여유를 확인합니다. 부족하면 계산을 시작하지 않고 CPU0·RAM4GiB만 유지하며 기다린 뒤 원래 예약을 다시 확보합니다. 이 보호는 이미 실행 중인 외부 벤치마크의 출력은 변경하지 않습니다.

### 하플로태깅 임시 조각 자동 정리

- 원인: 염색체별 병렬 WhatsHap BAM을 완성 BAM으로 합친 뒤에도 조각이 남아 같은 데이터를 중복 보관하던 문제를 수정함.
- 정리 시점: 완성 BAM·인덱스·태그 목록·단계 QC PASS를 확인한 후, 같은 입력의 완료 실행 및 더 오래된 실패 실행에서 재생성 가능한 임시 BAM만 정리함. 완성 BAM은 후속 전처리에 계속 사용하며 기존 최종 QC 정책을 유지함.
- 보존 조건: 열린 파일·검증 불일치·새 실행·symlink는 정리하지 않음. 원본·완성 BAM·QC·실행 장부·로그·태그 목록을 보존함. 정리 목록·바이트·QC 해시는 `haplotagged/parallel_haplotag/shard_cleanup.json`에 기록함.
- 실행: `06_haplotag.sh` 및 내부 Python 실행에 통합함. 별도 운영 스크립트 없이 기본 적용함. 문제 재현을 위해 조각이 필요한 경우 `MFM_HAPLOTAG_KEEP_SHARDS=1`로 보존 가능함. 기존 pysam·samtools를 사용하므로 새 환경 설치는 필요하지 않음.
- 확인된 회수량: 2026-10-29 기존12개 작업의 중복 조각1,727,412,989,952바이트를 제거함. 원본이나 모델 입력을 삭제한 용량이 아님.

- 같은 sample ID의 별도 CPU/GPU 비교: 공유 CPU/RAM 풀은 유지하고 실제 출력 경로로 작업 소유권을 구분함. 다른 출력의 GPU 비교를 기존 CPU 작업으로 오인하거나 샘플 RAM을 중복 귀속하지 않음. 같은 실제 출력은 실행번호가 달라도 중복 계산 방지 대상으로 보존함.
- 임대 반환 경합: CPU 예약 파일 반환을 `pool.lock` 안에서 수행함. 과거 manager의 반환과 겹쳐 사라진 다른 임대 파일은 스냅샷에서 건너뜀. 본인 예약 누락이나 장부 JSON 손상은 계속 오류로 처리하여 예약 없는 계산을 허용하지 않음.

### 2026-10-30 실제 CPU 활용 개선

- 최신 운영 RAM 상한1.65TB·예약1.485TB, 실제 CPU quota68코어를 유지합니다. 위의 임시 WGS01 한 개 제한은 아래 자동 산정으로 대체합니다.
- BAM 쓰기 동시 수는 현재 디스크 여유에서150GB를 남기고 작업당350GB로 산정합니다. 개별 쓰기 직전에는 원본 크기×2+50GB와 최소350GB도 확인합니다. 이미 진행 중인 쓰기는 종료하지 않으며 이후 신규 진입에 반영합니다. 선행 다운로드 동시1개·최소500GB 여유 조건은 유지합니다.
- 정렬된 BAM의 읽기 QC는 쓰기 슬롯과 분리하고 CPU/RAM 범위에서 진행합니다. 쓰기 대기 단계는 CPU0·RAM4GiB로 반환하며, 쓰기 직전 최신 정책으로 CPU/RAM·공간을 다시 확보합니다.
- 유휴 대여 기본값은 기존 최대 한 풀입니다. 현재 운영 정책은10초 실측 평균 유휴·nice19·실제 kernel quota 조건에서 낮은 우선순위 대여 최대136, 전체 논리 예약 최대204로 허용합니다. 이 수는 실행 CPU 수가 아니며 실제 전체 CPU 상한은68입니다. RAM 예약은 모두 합산합니다.
- 명시 정책 relend_fixed_prepare에서는 CPU/RAM이 고정된 WGS01의 계산·스레드수·RAM을 유지하고 해당 프로세스 트리와 기존 스레드를 nice19로 낮춘 뒤 정상 CPU 예약을 반환할 수 있습니다. 실제 유휴와 대여 상한을 확인하며, 후속 가변 자원 단계는 기존 정상 입장 검사를 사용합니다. 별도 프로젝트는 정지하지 않습니다.
- 실제 전처리 CPU와 공유 cgroup 전체 CPU를 구분해 측정합니다. I/O 및 직렬 구간, RAM·디스크 대기로 실제 사용량은 변합니다. 관련144검사 통과. 원본/전달의 전처리 내용과 cleanup 차이는 유지합니다.

### 2026-10-30 BAM 전수 집계와 입력 검증 병렬화

- BAM flag·HP·PS 전수 집계는 인덱스 없이 한 번 읽습니다. HTSlib의 BGZF 해독, 완전한 레코드 배치 생산자, 네이티브 집계 큐를 겹쳐 실행합니다. QC 실패·secondary 우선 supplementary·duplicate·미배치·HP/PS 존재 여부를 유지하며 중간 손상과 입력 변경은 실패로 처리합니다.
- I/O 보조 스레드를 포함한 모든 네이티브 스레드가 실제 CPU grant 안에 들어갑니다. 배치8MiB·레코드64MiB·대기 큐와20% 여유로 사전 RAM을 계산합니다. 8CPU에서 약1.218GB,68CPU에서 약7.459GB이며 기존 단계 예약과 전처리 RAM 상한1.65TB/예약1.485TB도 유지합니다.
- `MFM_RECORD_STATS_SCAN_MODE=auto|legacy|single_pass`이며 기본auto입니다. BAM은 새 경로, CRAM은 기존 samtools 경로입니다. C++17 컴파일러와 현재 Python 환경의 HTSlib 개발 파일이 없거나 실행 불가하면auto는 기존 경로를 씁니다. 최초 빌드는 사용자cache의 `mfm_native`에 보관하며 `MFM_NATIVE_CACHE_ROOT`로 변경할 수 있습니다. 입력 손상 오류는 정상 결과로 대체하지 않습니다.64MiB보다 큰 레코드는auto에서 기존 경로로 처리합니다.
- 실제 HG02602 HiFi6만 레코드796,605,229바이트에서 같은8CPU 예산·순서를 번갈아3회 비교했습니다. 기존 세 번 집계 중앙값4.723초→한 번 집계1.681초(2.809배), 새 방식의1CPU3.673초 대비2.185배입니다. 카운트는 모두 일치했습니다. 새 방식 RSS 약50.9MB·실제 스레드8개입니다. 이 수는 해당 집계 단계의 측정이며 전체 전처리 배율이 아닙니다.
- 명시 정책 `cpu_backfill.validation_phases`에서는 WGS02 입력 검증도 실제 유휴 CPU를 nice19로 빌릴 수 있습니다. 하위 병렬 검증이 CPU를 반환한 뒤 다음 묶음에서 최초 grant·대여 전체 상한·전체 논리 상한·RAM 예산 안에서 다시 확보합니다. 재정렬/compact는 대상에 포함하지 않습니다.
- 이미 import한 실행 중 프로세스를 종료하지 않습니다. 새 프로세스부터 새 집계를 사용하며, 대기 중인 입력 검증은 새 자원 정책으로 진입합니다. 원본과 전달 코드의 cleanup 차이는 유지합니다.

#### public VCF의 염색체 병렬화 추가

public VCF03은 염색체마다 기존 추출·정규화 함수를 독립 임시 경로에서 실행한 뒤 concat·전역 정렬·index를 수행합니다. 염색체 경계에서 나누며 GT·PS·중복·INDEL·sample alias를 유지합니다. 실제 bcftools의 직렬/병렬 전 필드·통계 일치와 REF 오류 시 부분 출력 거부를 검사했습니다. `MFM_PUBLIC_VCF_JOBS=auto|양수`이며1은 기존 직렬 경로입니다.

CPU는 worker당 native pipe2개와 부모1개, RAM은 worker당2GiB·부모2GiB에20%여유를 계산합니다. 정렬기 기본 버퍼768M와 FASTA/BCF 버퍼를 포함하며 실제 lease·Slurm 예산 안에서 동시 수를 정합니다.27CPU이면13염색체·약36.08GB의 피크 예측입니다. 기존 예약 이하일 때만 피크 예측으로 RAM 예약을 낮추고 초기 예약과 산정 근거를 lease에 보관합니다. RAM 증액과 CPU 증설은 하지 않습니다.

명시 유휴 대여 정책은 WGS02 입력 검증과 WGS03 VCF 준비의 진입을 지원합니다.03은 단계 중 CPU를 늘리지 않으며02만 하위 묶음 사이의 대여 CPU 재확보를 사용합니다. 실제 사용량과 논리 예약을 구분하고 현재 전체 kernel CPU68·전처리 RAM1.65TB/예약1.485TB를 유지합니다.


### 2026-10-30 public VCF 준비와 haplotag 예약 분리

- 적용 범위: public VCF와 제공된 위상을 사용하는 HiFi long-read 경로임. 실행은 기존 `run_pipeline.sh run`을 사용하며 별도 제어 프로그램·새 설정 파일은 필요하지 않음.
- 입장 RAM: WGS03 VCF 준비와 WGS06 haplotag의 기존 공통128GiB 예약을 분리함. VCF 준비는 `ceil((부모2GiB + 동시 worker×2GiB)×1.2)`로 산정하며, 동시 worker 상한은 `min(22, max(1,(CPU-1)//2))`임. 68CPU·상염색체22개의 최대55.2GiB이며, 실제 grant와 염색체 수가 작으면 실행 시 추가 감축함. haplotag128GiB는 유지함.
- 완료 반환: 염색체별 변환 worker와 정렬·index가 모두 종료된 후 태그 생성·스트리밍 검증 전에 CPU최대2개·RAM4GiB로 감축함. 완료본 재사용 경로에도 동일 적용함. 다음 haplotag는 내장 단계 큐에서 별도로 자원을 확보함.
- 호환성: 과거 `MFM_LONGREAD_MEMORY_PLAN`의 `public_vcf_and_haplotag` 계획도 WGS03에서 새 산식으로 해석함. 새 계획은 `public_vcf_prepare`와 `haplotag`로 구분함. 입장 후 예산 초과 증액 거부·본인 예약 누락 오류·전체 RAM 보호는 유지함.
- 검증: 원본/전달용 각각55개 관련 검사 통과함. 실제 bcftools의 직렬/병렬 seq/tag 레코드·GT·PS 동일성, 태그 생성 전 자원 반환, 캐시 재사용, 구형 계획의 작은 RAM 입장, CPU/RAM 상한을 확인함. 전체 대용량 샘플의 처리시간 개선 폭은 별도 실측이 필요함.

#### haplotag 병렬 실행의 유휴 CPU 진입

명시 정책 `cpu_backfill.haplotag_phases`에서 WGS06도 실제 CPU·RAM 여유가 확인되면 nice19로 진입합니다. 기존 염색체별 WhatsHap 병렬화와 태그·미배치·출력 순서를 유지합니다. 대여 단계는 최초 grant를 넘어 증설하지 않으며 남은 염색체에 따라 CPU를 반환합니다. 단계 RAM128GiB 예약은 종료까지 유지합니다. 실제 WhatsHap의 직렬/병렬 BAM·CRAM 동등성과 기존 부분 결과 보존을 재검증했습니다.

### 2026-10-30 전처리 우선 배분 보완

운영 공유 정책의 `cpu_backfill.nice=0`은 새 전처리 대여를 일반 우선순위로 실행합니다. 기본값19는 유지하며 기존 nice19 manager는 권한 없이 상향하지 않고 실제값을 기록합니다. 운영에서 `relend_fixed_prepare=false`로 두어 정상 전처리의 우선순위를 새로 낮추지 않습니다. 실제 cgroup CPU68·RAM1.65TB/예약1.485TB 및 디스크 보호는 유지합니다.

명시한 HiFi RAM 모델에 동적 재정렬과4GiB 준비 입장 예약이 선언돼 있으면 WGS10 준비 최소 CPU만2로 시작합니다. 실제 정렬은 별도 phase에서 기존 CPU/RAM과 최소 정렬 스레드를 다시 확보합니다. legacy·동적 phase 비활성·일반8CPU 계약은 유지합니다. 준비 대기 중 계산 자식이 없는 두 manager만 재진입시켰으며 계산 중인 작업은 유지했습니다.

사용자가 GLM 수정 금지를 명시했으므로 별도 GLM의 CPU·우선순위·affinity·설정은 변경하지 않습니다. 공유 CPU에서 전처리만68코어를 독점한다고 보고하지 않습니다.


### 2026-10-30 haplotag 후반 자원 분리와 GPU 검토

- 적용: 단일 샘플·내장 WGS06 공유 lease 경로에 적용함. 기존 run_pipeline.sh를 사용하며 추가 운영 제어기는 필요하지 않음. 여러 샘플이 하나의 lease를 공유하는 외부 호출에는 후반 예약을 적용하지 않음.
- RAM: 염색체별 WhatsHap 계산의128GiB는 유지함. 모든 worker가 종료한 뒤 unmapped 복사·BAM 결합·태그목록 스트리밍 처리에는4GiB만 남김. 일부 살아 있는 염색체의 RSS만으로 전체 계산 피크를 낮추지 않음.
- 인덱싱: samtools index 직전에 설정된 output_threads+1 범위에서 공유 CPU를 재확보함. 실제 예약 CPU가 n이면 index 보조 thread는 최대 n-1임. RAM은 부모4GiB+예약CPU당64MiB로 산정함. 태그 검사에는 최대2CPU, 태그 집계에는1CPU를 예약하고 RAM4GiB로 반환함. 각 전환은 HAPLOTAG-PHASE 로그와 lease의 haplotag_resource_events에 기록함.
- 유휴 대여: 명시 haplotag_phases 정책에서 최초 대여 grant·낮은 우선순위·대여 전체 상한·논리 예약 상한·공유 RAM 안에서만 후반 CPU를 재확보함. RAM이 부족하면 계산 전에 대기하며 본인 lease 누락은 오류로 처리함. 전체 CPU quota68·RAM1.65TB/예약1.485TB는 유지함.
- HiFi 재정렬: adaptive 계획의 준비 진입 최소CPU2를 original stage10·전달 layout stage08에 동일 적용함. 이후 실제 집계·FASTA 준비·정렬은 기존 RealignPhases에서 전체 풀을 활용해 별도 확보함. 정렬의 최대CPU를2로 제한하는 설정이 아님. adaptive 비활성·명시 초기CPU·short-read 설정은 보존함.
- GPU 확인: WhatsHap 공식 문서·Parabricks4.8 공식 도구 목록·설치 바이너리에서 현재 HP/PS/PC 및 read 보존 계약을 그대로 대체하는 GPU haplotag 구현은 확인하지 못함. Clair3 GPU의 내부 in-memory haplotag는 별도 변이 검출 흐름이며 논문에서도 WhatsHap과 태깅 결과 차이가 보고됨. 따라서 현재 WhatsHap을 유지함.
- 근거: [WhatsHap guide](https://whatshap.readthedocs.io/en/latest/guide.html), [Parabricks4.8 도구 목록](https://docs.nvidia.com/clara/parabricks/about-parabricks), [Clair3 GPU 논문](https://pmc.ncbi.nlm.nih.gov/articles/PMC13141153/).


## CPU·GPU HAP 병렬 및 입력 준비 개선 (2026-10-30)

- CPU HAP 동시수는 샘플 수의 사전 분할이 아닌 실제 단계 CPU/RAM으로 결정한다. 자동 정렬 thread·sort 보조 thread를 조합하며 명시 HAP/스레드 설정은 유지한다. 원본 로컬 실행기도 기본 병렬 허용이며 명시 `MFM_HAP_PARALLEL_ALLOWED=0`은 보존한다.
- 단일 샘플 Illumina 준비는 4GiB로 진입한 뒤 `input_prepare`에서 CPU와 `4GiB+2GiB/CPU`를 재예약한다. 실제 입력 통계 모듈은 batch·QNAME·decode·worker 피크를 다시 검증한다. 단계 예약 비활성/다중 샘플 내부 jobs의 기존 진입 예산은 유지한다.
- GPU 입력 준비·대기·계산·출력은 HAP별 자원 상태를 합산한다. 장치 대기는 CPU1/RAM4GiB, 공간 대기는 CPU0/RAM4GiB이며 GPU와 host CPU/RAM·쓰기 예상량을 함께 확보하지 못하면 장치를 즉시 반환한다. 도구 계산 host 기본 128GiB와 10% 여유는 실제 계산 동안 유지한다.
- 같은 `all_to_both` 입력은 잠금·입력 지문·완료 장부가 맞는 FASTQ만 HAP 사이에 공유한다. 부분 생성/손상/입력 변경은 재사용하지 않는다. 4MiB 배치, FASTQ 행 64MiB 상한과 실제 grant 안의 pigz worker로 준비한다. 압축 worker가 부족하면 기존 gzip 경로를 사용한다. RAM 예약 4GiB에는 최악 레코드·batch·producer·압축 버퍼의 여유를 포함한다.
- GPU BAM은 RG·OD duplicate를 한번 복원한 비압축 BAM으로 CRAM writer에 전달한다. 이미 생성된 구형 paired shell에는 기존 SAM 출력을 유지한다. RNA coverage와 깊이 QC는 GPU 배열을 재사용하며 Zarr Blosc/Zstd·독립 CPU source trace 계약은 유지한다.
- 임시 쓰기는 입력 준비 `2×입력 BAM+8GiB`, GPU 출력 `3×압축 FASTQ+32GiB`를 보수 예상하고 디스크 여유·다른 HAP의 남은 예약·정책 여유를 대조한다. 실제 quota/cgroup I/O의 강제 격리와 다르며 도구 오류는 계속 실패 처리한다.
- 현재 실행 중인 Python/셸을 재시작하지 않는다. 이미 확정된 단일 HAP 루프는 유지하고, 아직 생성하지 않은 다음 HAP/다음 단계·새 실행에서 새 모듈을 읽는다. 두 HAP 병렬 계획은 새 재정렬 진입부터 적용한다.
- 실제 HiFi 2,048개/87.64MB 부분 입력 시험: 공유 풀 12CPU 조건에서 변환·압축·쓰기 중앙값 0.706→0.087초(8.08배), 전수 해독 SHA256 일치. 이 수치를 전장 전처리 배율로 해석하지 않는다. GPU 점유 중인 전장 정렬의 새 실측은 아직 없다.


## GPU 모드의 Parabricks 도구 선택 (2026-10-30)

- 실행: 기존 `run_pipeline.sh run --selection FILE --config FILE --executor local`에 `--acceleration-mode gpu --gpu-devices 0,1,2,3 --gpu-wait-seconds 86400`을 추가하거나 config의 `acceleration.mode`를 `gpu`로 지정함. 별도 샘플 제어 Python은 필요하지 않음. GPU 실행환경은 기존 설치 안내를 따름.
- 개별 선택: `--bam-backend cpu|parabricks` 또는 config의 `bam_backend`로 BAM 도구만 선택 가능함. 명시한 개별 backend는 모드 기본값보다 우선함. 기존 config에 CPU backend가 명시되어 있으면 GPU 모드에서도 해당 선택을 유지하므로 함께 확인함.
- RNA GPU: WGS-only 선택은 사용하지 않는 RNA 호환성 검사를 요구하지 않음. RNA pair는 기존 strict/experimental 규칙을 유지함. 현재 프로필과 다른 GPU STAR 정책을 시험할 때만 `--rna-gpu-policy experimental`을 명시하고 CPU 대조를 수행함.

| 과정 | GPU 모드의 실제 구현 |
|---|---|
| WGS 재정렬 | Parabricks `fq2bam`/`minimap2` 사용함. supplied phasedVCF·원본 duplicate 정책을 유지하고 추가 BQSR·변이검출을 수행하지 않음. |
| RNA BAM→FASTQ | Parabricks `bam2fq` 사용함. Parabricks 제공 CPU 도구임. paired suffix만 복원하고 orphan/unpaired/QC-fail 처리 및 read 이름·서열·품질 검증을 유지함. |
| RNA 정렬 | 기존 Parabricks `rna_fq2bam` 경로 사용함. |
| RNA 중복 표시 | GPU STAR와 GPU BAM backend를 함께 선택하면 `rna_fq2bam` 통합 중복 표시를 기본 사용함. 외부 queryname 재정렬·별도 markdup 호출을 줄임. coordinate baseline·optical distance100·duplicate 보존·metrics와 독립 flag QC/index를 유지함. `RNA_GPU_MARKDUP_MODE=separate`는 기존 분리 경로임. |
| depth·QC·서열 packing | 기존 CuPy 선택을 사용함. 실제 1bp signal·uint32·ACGTNP4bit·Zarr Blosc/Zstd5 계약을 유지함. |
| VCF consensus·haplotag·입력 무결성 | 검증된 CPU 경로 유지함. Parabricks aggregate bammetrics는 기존 독립 source trace 및 1bp track을 대체하지 않음. |
| BAM 압축·CRAM CRC | 별도 `--io-backend nvcomp` 선택임. 작은 실제 시험에서 CPU보다 느려 GPU 모드의 자동 기본값은 CPU임. 전장 비교 결과와 구분함. |

- 중복 CPU 정렬 제거: GPU가 coordinate 정렬한 WGS BAM은 header를 확인하고 RG/OD만 복원한 뒤 바로 CRAM으로 인코딩함. 출력 정렬·index·독립 QC 검증은 유지함. 이미 실행 중인 구형 shell은 기존 경로로 완료함.
- 자원·대기: 도구별 thread/host RSS는 실제 CPU·RAM 예약 안으로 제한함. GPU0–3 공유 풀은1초마다 장치를 확인함. 외부 사용 중인 GPU를 신규 점유하지 않으며 CPU/RAM·디스크 조건도 충족해야 계산을 시작함. 외부 프로젝트가 나중에 같은 GPU를 재사용하는 것까지 강제 차단하지는 못함.
- GPU writer 헤더: 중복 `@HD`가 생성되는 실제 사례를 확인함. coordinate 헤더 존재를 확인하고 중복 HD만 `samtools reheader -P`로 제거함. read 본문·HP/PS·순서를 그대로 보존하는 회귀 검사를 적용함. 중복 헤더가 없으면 추가 복사를 수행하지 않음.
- 실물 검증: CPU 제공 `bam2fq`는 samtools와 read 이름·서열·품질·제외 수 대조를 통과함. `tests/test_parabricks_bam.py`에 실제 GPU `bamsort/markdup`와 Picard의 동일 read·duplicate flag·HP/PS 대조가 있음. `MFM_RUN_PARABRICKS_GPU_TESTS=1`로 별도 실행함. GPU가 점유 중이면 검증도 장치를 기다림.
- 운영 검증 표식: 현재 세 비교 샘플은 `MFM_PARABRICKS_MARKDUP_VALIDATION`의 실제 GPU 통과 판정·backend SHA256을 확인한 뒤 RNA 중복 표시를 수행함. 진행 중에는 기다리고 실패·코드 변경 시 출력 생성을 차단함. 전장 CPU/GPU 정확도와 속도 비교는 전체 샘플 완료 후 판단함.
- HG00126 비교 재개:128GiB host RSS guard 초과를 확인해 해당 비교 설정만192GiB로 조정함. 기존 HAP당 약220GiB 예약 모델과 공유 RAM 예산을 적용하며 `run_pipeline.sh resume --from-stage 10` 후 내장 run으로 재개함.
- NVBIO 검토: CUDA/C++ 연산 라이브러리로 SW/Gotoh·FM-index·nvBowtie short-read 등의 구성요소를 제공함. 현재 HP/PS/PC·주어진 phasedVCF 계약의 완성된 haplotag 대체는 확인하지 못함. 새 알고리즘 구현과 현재 CUDA 호환성 검증이 필요하므로 운영 backend에는 추가하지 않음.
- 공식 근거: [bam2fq](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/bam2fq), [bamsort](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/bamsort), [markdup](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/markdup), [NVBIO](https://nvlabs.github.io/nvbio/), [NVBIO 소스](https://github.com/NVlabs/nvbio).


### 정렬 하나에 GPU3장과 별도 로컬 임시 공간 사용 (2026-10-30)

- 설정: `pipeline_env`의 `MFM_GPU_ALIGNMENT_GPUS=3`, `MFM_GPU_DEVICES=0,1,2`, `CUDA_VISIBLE_DEVICES=0,1,2`를 사용함. 현재 비교 실행은 `MFM_GPU_HAP_JOBS=1`·`WGS_HAP_JOBS=1`로 HAP를 순차 처리함. 같은3장 풀에서 정렬 하나만 실행 가능하며 준비된 다음 정렬은 GPU 반환 후 자동 투입함. 샘플의 모든 CPU/RNA 단계를 통째로 직렬화하는 규칙과 구분함.
- 예약: 기존 공유 GPU lock 안에서3장을 모두 확보하거나 아무 장치도 확보하지 않음. 외부 사용 장치·Slurm/CUDA 공개 범위는 보존함. Parabricks `--num-gpus`·native/container 공개 UUID 개수를 맞춤. 기본값1장과 기존 단일GPU coverage/QC/중복표시는 유지함.
- 메모리: 여러 GPU의 VRAM을 단일 메모리로 합치는 기능이 아님. Parabricks가 작업을 여러 GPU로 분산함. 현재 세 비교 설정은 host 한도384GiB와 입력·버퍼·10%여유를 반영한 HAP당432GiB 예약을 사용하며 전체 CPU68·RAM1.65TB/예약1.485TB 및 실제 할당을 우선함.
- 임시 저장: `MFM_GPU_SCRATCH_ROOT`를 지정하면 새 FASTQ 캐시·GPU BAM·중간 CRAM을 그 파일시스템에 생성함. `/tmp`가 tmpfs인 서버에서는 RAM 파일시스템을 대형 임시 공간으로 선택하지 않음. 현재 비교 실행은 `/home/aibio_1/.cache/mfm_gpu_261030`의 별도 overlay 공간을 사용함.
- 경로 보존: 원래 출력 위치에는 새 로컬 중간 CRAM을 가리키는 링크만 생성함. 기존 파일은 이동하거나 덮어쓰지 않음. 모델 입력 Zarr·QC·로그·기존 CPU 완료본 경로는 유지함. 중간 CRAM 링크는 로컬 파일에 의존하므로 전장 비교가 끝날 때까지 보존함. 다른 노드에서 재개할 때에는 해당 경로 접근성을 별도로 확보해야 함.
- 검증: 원본/전달 각77개 검사 중72개 통과·GPU 수치5개는 별도 opt-in으로 생략함. Docker 다중장치 인자의 CSV 해석을 포함한 관련28개 검사도 각각 통과함. 생성된 실제 shell을 두 HAP에 실행해 로컬 CRAM 링크·read·HP/PS·duplicate·입력 공유를 검증함. 실제 HG00126 개인 HAP 구간 유래64개 입력의1GPU/3GPU 결과111개 alignment에서 전체 SAM 필드·태그가 같았음. 작은 입력의 시간은 전장 속도 배율로 해석하지 않음.
- CPU 대기 계산: 이전 하위 worker가 종료한 `realign_phase=:waiting` 구간은 정상 CPU 예약을0으로 계산함. 대기 부모의1CPU씩 때문에 GPU 준비 진입이 막히던 문제를 보완함. 실제 계산 예약과 전체68CPU 한도는 유지함. 현재 비교 설정은 GPU host CPU 상한68·최소 align thread1로 가용CPU에 맞춰 실행함.
- GPU 유휴 대여: 공유 정책의 `cpu_backfill.gpu_realign_phases=true`에서만 GPU 경계의 대여를 허용함. 기존10초 실측·nice·전체논리204/대여136·실제CPU68·RAM 한도를 검증함. GPU 대기 후 재확보도 최초 대여 CPU 상한을 넘지 않음.
- 준비 BAM과 잠금: Illumina route/collate/paired BAM 준비도 GPU scratch에 새로 생성함. CPU 기본 경로는 유지함. 공유 FS의 blocking flock에서 일시 EAGAIN이 발생하면 재시도하며 권한 등 다른 오류는 보고함. 비교용 원본과 완료 결과는 보존함.
- 남은 저장 공간: 로컬 임시 저장은 Lustre 쓰기 대기를 줄이지만 최종 모델 입력을 저장할 Lustre 공간까지 늘리지는 않음. 운영 데이터 삭제나 타 프로젝트 제어는 수행하지 않음. Slurm에서는 요청 GPU 수 이상의 `slurm.gpus_per_task`와 노드별 scratch 접근성을 설정해야 함.

### 2026-10-30 재정렬 RAM 활용 보완

공유 정책의 `cpu_backfill.realign_phases=true`는 재정렬의 완료 경계에서 실제10초 평균유휴 CPU를 대여합니다. 고정 단계와10초 안정화 시간을 공유하며 전체논리204·대여136·실제cgroup68코어 및 RAM예약1.485TB를 함께 지킵니다. 실행 중 계산의 스레드와 우선순위는 바꾸지 않습니다. 내부 자원대기는 감시기의 `internal_waiting_resources`에 따로 표시됩니다.

`realign_sort_memory_bytes`는 자동 sort 버퍼의 최대치입니다. 현재64GiB/스레드에서 예산별로 절반씩 줄여 선택하고 실제 HAP/스레드별 피크를 시작 전에 예약합니다. 큰 버퍼로 HAP동시수가 줄면 작은 버퍼를 선택하며 명시된 더 큰 기존값은 보존합니다. RAM운영상한1.65TB·10%여유 및 디스크 안전대기는 유지합니다. GLM 코드·프로세스·자원설정은 변경하지 않습니다.

재정렬의4GiB동적준비는동적계획을검증한경우에만2CPU로유휴입장하며실제계산피크는다시예약합니다. RAM피크없는계산입장을허용하지않습니다. 큐의장부읽기는원자JSON스냅샷을사용하고실제할당기만pool.lock에서CPU/RAM을확보합니다. 감시기는같은lock안에서호스트/PID/시작시간소유권을검증하며, 별도보조잠금파일의고아경합으로조용히종료하지않습니다.


### GPU 전장 비교의 예약 예외 및 4장 실행 (2026-10-30 추가 승인)

- 선택 설정: `run_id`가 `gpu_compare_`로 시작하고 Parabricks backend를 사용하는 비교에만 `pipeline_env.MFM_GPU_COMPARE_UNRESERVED=1`을 허용함. CPU·RAM 입장 검사·대여 대기·단계별 예약을 생략함. 일반 CPU/배치 실행의 정책은 유지함.
- 실제 한도: 스레드 상한은 affinity·cgroup·Slurm 중 작은 CPU 한도에 맞춤. 현재68코어이며 고정 할당 또는 실사용68코어를 뜻하지 않음. 실제 host RSS 보호·커널 메모리 상한·GPU 독점 잠금·full QC·source trace는 유지함. 비교 도구 host RSS 상한384GiB를 적용하며 RAM을 선점 예약하는 값과 구분함.
- 네 장 사용: 현재 세 비교 설정은 `MFM_GPU_DEVICES=0,1,2,3`, `CUDA_VISIBLE_DEVICES=0,1,2,3`, `MFM_GPU_ALIGNMENT_GPUS=4`를 사용함. HAP 정렬은1개씩이며 Parabricks가4개 GPU에 작업을 분산함. 이전3장 설정에 우선하는 최신 비교 설정임.
- 실행 순서: CPU 준비·출력·QC worker는 동시 진입함. GPU 도구 단계에서 공용 풀의4장 독점 잠금으로 정렬 하나씩 실행함. 준비된 다음 GPU 작업은 반환 후 진입하며 CPU 출력 때문에 전체 샘플을 직렬화하지 않음. `run_pipeline.sh resume/run`으로 재개함.
- 재사용: 완료 증명이 있는 FASTQ를 재사용함. NA18959 82.5GB·HG00126 105.37GB 완성 입력을 보존함. 기존 CPU 완료본·raw·중간자료는 전장 비교가 끝날 때까지 보존함.
- 시간 기록: `resource_mode=unreserved_gpu_comparison`의 CPU audit은 예약CPU/RAM0·실제CPU초·벽시계시간·RSS를 기록함. 단계 이벤트의 계산 계획값은 `reservation_applied=false`와 함께 기록되며 실제 예약량으로 해석하지 않음. 샘플 잠금 대기·중단된 시도·완성 입력 재사용은 정렬/전체 처리 시간과 구분해 비교함.
- 종료 처리: Runtime의 SIGTERM/INT는 분리된 Parabricks 프로세스 그룹까지 전달하고 종료를 기다림. 네 장 전환 과정에서 이전3장 자식이 남던 경로를 실제로 확인해 보완했으며 별도 자식 종료 검사로 확인함.
- NVIDIA 입력 경로 검증: 설치된 `nvCOMP nvlzcat`의 GPU Gzip 압축 후 CPU 해독 바이트 일치를 확인함. 압축 부분의 GPU 사용 가능성을 뜻하며 전체 BAM→FASTQ GPU 변환이나 전장 속도 향상을 입증한 결과는 아님. `minimap2 --in-bam` 소량 실행은111개 alignment를 생성했지만 HP/PS 태그가 소실돼 현재 FASTQ 경로를 대체하지 않음.


### GPU HiFi의 전 배치 recovery 보완 (2026-10-30)

- 확인 현상: HG00126 전장 HAP1의8,809,923read/chunk1000에서GPU seeding/chaining recovery8,810회를 확인함. GPU가 메모리를 점유해도 실제 seed/chain 계산이 CPU 복구 경로로 진행될 수 있어 GPU 개수·VRAM과 실제 연산을 구분함.
- GPU HiFi 기본값: `map-hifi`의 완료 FASTQ 앞1,000개 길이를 읽어 p95를 산출함. `--estimated-read-length=max(15000,p95)`, `--estimated-anchor-per-read=max(7000,2×추정길이)`로 NVIDIA 공개 GPU 버퍼 설정을 보완함. read·품질·HP/PS·VCF·window 필터를 변경하지 않음. Illumina/ONT 자동 보정은 추가하지 않음.
- 명시 설정: `MFM_GPU_MINIMAP_ESTIMATED_READ_LENGTH`와 `MFM_GPU_MINIMAP_ESTIMATED_ANCHOR_PER_READ`를 지정하면 우선함. `MFM_GPU_HIFI_AUTO_PROFILE=0`은 자동 산정을 끔. `MFM_GPU_MINIMAP_LARGE_CHAINING_SIZE`는 사용자 명시 시에만 적용함. 장문 chaining 기준의 기존 기본값은 유지함. `MFM_GPU_MINIMAP_VERBOSE=1`로 GPU 커널 시간 등 진단 로그를 남김.
- 입력 보존: 기존100GB 이상 FASTQ를 재생성하지 않고 앞부분의 길이 통계만 읽음. `gpu_input_performance_profile`에 표본수·평균·p95·최대길이·실효 버퍼값을 저장함. 앞부분 표본이 전장 anchor 분포를 보장하지 않으므로 나머지 입력의 recovery도 계속 기록함.
- 실제 검증: 같은 HG00126 전장 참조·실제1,000read에서 버퍼만 변경함. 기본 recovery1회·GPU seed/chain0초 → 보정 recovery0회·GPU seed0.02초/chain0.07초이며4,594개 alignment의 전체 SAM 필드·태그가 일치함. 장문 cutoff를 유지한 분리 검증도 통과함. 작은 입력의 초기화 포함 시간은79.4→94.8초로 증가했으므로 전장 속도 향상을 보장하는 결과로 해석하지 않음.
- 실제 자원: 현재 비교 실행의 host RSS 보호384GiB와 전체RAM1.65TB/실제cgroup한도를 유지함. 버퍼 보정으로 host/GPU 메모리가 증가하므로 실측함. GPU 비교의 CPU/RAM 예약 대기는 생략하되 실제 RSS 보호는 유지함.
- 진단 메타데이터: 로그에 있으면 `gpu_seeding_chaining_recovery_count`, `gpu_seed_kernel_seconds`, `gpu_chain_kernel_seconds`를 기록함. recovery0만으로 전체 전처리가 GPU에서 수행됐다고 해석하지 않음.
- 동시 실행: 전체 샘플 잠금을 제거하여 CPU 출력/QC 중 다른 준비 완료 샘플이 GPU에 진입함. 기존4GPU 독점 잠금과 실제 CPU quota68을 유지함. 표식 JSON은 기존 고유 tempfile 원자 쓰기를 사용함. CPU 배분 검사는 고정8코어 가정 대신 선언된 최소CPU 이상과 동시 CPU/RAM 총합을 검증함.

### GPU compact 병렬화·자원 반환

- 공유 GPU 계산은 장치별 예약 비율·worker 수를 기준으로 분산함. 허용 GPU 목록·독점 정렬 예약·다른 compute PID 보호 유지함. 16MiB 이하의 드라이버 잔여 메모리는 실제 compute PID 없음·활용률0을 확인한 경우만 유휴로 취급함.
- GPU compact는 spawn worker에서 염색체 한 개 처리 후 worker를 교체함. 완료된 worker의 CUDA context와 예약을 반환함. CPU 기본 경로의 worker 재사용은 유지함.
- `MFM_COMPACT_ALIGNMENT_THREADS`로 CRAM 해독 thread 지정 가능함. CPU/GPU 모두 실제 CPU 예산에서 자동 계산함. HAP별 CRAM은 순차 처리하므로 동시 handle1개이며, 68CPU/22worker에서는 계산1+해독보조2개로 총66CPU 이내임. 추가 해독 버퍼64MiB/worker를 RAM 예측에 포함함. Blosc thread와 명시 해독 thread도 같은 CPU 예산에 포함하며 실제 사용률·속도 향상은 별도 실측함.
- source trace는 provenance와 실제 파일이 같음을 확인한 뒤 원래 CRAM 링크 경로를 사용함. 링크 옆 `.crai`를 유지하며 다른 파일로 바뀐 링크는 거부함.
- 실제 NA18959 4개 염색체·각 HAP 2개1Mb창에서232개 배열 동일·GPU4장 분산·worker 종료 확인함. 작은 구간 전체시간37.21→38.02초로 전체 속도 향상은 미확인임. GPU 메모리 점유·활용률·CPU/I/O 구간을 구분함.


### FASTQ 동시 공급·태그 복원·HAP별 Zarr 병렬 처리

- GPU 입력 : `MFM_GPU_INPUT_MODE=auto` 기본임. 생산자와 정렬의 CPU 예산이 함께 확보되면 FIFO로 FASTQ를 공급하면서 기존 gz 캐시를 생성함. HAP2는 검증 완료 캐시를 재사용함. `stream`/`cache` 명시 선택 가능함. CPU가 작으면 auto는 기존 cache 경로를 사용함.
- 입력 계약 : read·mate 순서·SEQ·QUAL·HP/PS/OD/XU 보존함. 첫 profile read도 다시 공급함. paired FIFO는 독립 bounded writer를 사용함. source/소비자 실패를 전파하고 불완전 캐시에 완료 장부를 쓰지 않음. FIFO/캐시 생성은 겹친 시간이므로 단계 wall에 다시 더하지 않음.
- 태그 복원 : `MFM_GPU_TAG_RESTORE=auto` 기본임. 기존 native 빌드 캐시의 HTSlib 루프로 RG/OD만 복원함. 컴파일 불가 시 기존 Python 경로 사용함. `native`는 빌드 실패를 오류로 보고하며 `python`으로 명시 비교 가능함.
- 비교 출력 : 현재 CRAM 출력은 아래의 염색체 병렬 경로를 사용함. 작은 예산의 직렬 경로에서는 CPU68 안에 복원 해독8·CRAM 인코딩 보조57·주프로세스를 배분함. CRAM 버퍼로 RAM을 산정함. CPU 대조군도 `gpu_compare_` 실행명과 `MFM_CPU_COMPARISON_CONTROL=1` 명시 시 같은 예약 예외·실제 CPU/RSS 보호를 사용함.
- GPU context : CuPy import 전에 공유풀에서 배정된 UUID 한 장만 공개함. 다른 카드에 숨은 CUDA context를 만들지 않으며 worker 종료 시 lease를 반환함.
- CPU 전달 : 정상 예약과 명시 비교가 같은 단계 환경 함수를 사용함. WGS03 염색체 병렬 수와 WGS06 WhatsHap CPU 예산이 실행 경로에 따라 누락되지 않도록 함.
- Zarr : 개인 재정렬 경로의 `MFM_COMPACT_HAP_JOBS=auto` 기본임. HAP별 독립 spawn 출력 후 부모가 기존 형식으로 결합함. 68CPU·Blosc2 기준17 contig×2HAP×2CPU=68CPU 예산임. 부모2GiB와 HAP별 기존 worker RAM 상한·해독 버퍼를 포함함. 실행 중 계산은 유지하고 새 작업에 반영함.
- 검증 : 실제 Parabricks single128/paired256 alignment의 파일/FIFO 입력 SAM 필드·태그·순서 동일함. 같은 실제16개 window의232개 배열에서 HAP 순차/병렬 결과 동일함. 작은16CPU 시험33.40→18.74초이며 전장 가속률과 구분함. CPU/GPU 정렬기 간 완전일치는 별도 전수대조 대상임.

- Illumina 성능 비교 옵션 : `MFM_GPU_CIGAR_ON_GPU=1`은 CIGAR GPU 계산, `MFM_GPU_BWA_PARTITION=1`은 2GPU씩 partition 실행임. GPU4장·충분한 보조CPU·두 옵션 모두 미지정이면 전량 대조를 통과한 동시 적용을 기본으로 사용함. 다른 장치 수는 기존 기본을 유지하며 0/1 명시값을 우선함. CIGAR만 명시하면 partition은 자동 활성화하지 않음. GPU0~3/4장에서는 partition 두 개를 사용함.


### GPU 정렬·CPU 출력 중첩 및 전수 trace 배치 처리

- HAP 연결 : GPU 기본 상한은 HAP 파이프라인2개임. HAP1의 CPU CRAM·인덱스·QC 중 HAP2의 GPU 정렬 진입을 허용함. GPU 계산은 기존 장치별 독점 풀에서 직렬화하며, HAP마다 GPU를 고정하지 않음. `MFM_GPU_HAP_JOBS=1`, `WGS_HAP_JOBS=1` 명시는 순차 실행을 유지함.
- CPU 인계 : GPU 반환 전에 CPU 출력 몫을 기존 단계 장부에서 먼저 교환함. 다음 HAP가 GPU를 기다리면 출력은 최대16CPU를 사용하며 나머지를 정렬 보조에 배분함. 68CPU에서는 출력16·다음 정렬 전체52가 가능한 예임. 실제 계산 코어 수와 구분함. 외부 CPU 대조군 완료를 기다리지 않음.
- 비교 장부 : 예약 예외 실행은 샘플 내부의 로컬 HAP 장부로 스레드 합계를 제어함. 공유 CPU/RAM 예약을 새로 걸지 않으며 기존 RSS 감시·cgroup 한도를 유지함. 정상 실행은 기존 공유 lease를 사용함. 실패/완료 HAP의 자원은 finally에서 반환함.
- CRAM 출력 : GPU 정렬 BAM을 reference별로 읽어 C++ RG/OD 복원과 CRAM 압축을 병렬 수행함. 원래 SQ 순서·동점 read 순서·태그·placed-unmapped 및 마지막 unplaced를 보존함. shard별 record 수를 BAM index 통계와 대조하고 quickcheck 뒤 결합함. CPU가 작거나 네이티브 사용 불가 시 기존 직렬 인코딩으로 처리함. 오류를 정상 fallback으로 숨기지 않음.
- 출력 자원 : fetch·태그 복원·인코더 main·codec·I/O를 worker당5CPU, 부모2CPU로 계산함. worker 수는 `min(구간 수, (CPU-2)//5)`임. RAM은 `(부모4GiB+worker당3GiB)×1.1`로 산정함. 성공한 자체 shard만 정리하며 실패 기록과 파일은 보존함.
- CPU HiFi/ONT : auto sort의 2스레드 상한을 제거함. HAP당 가용 CPU에서 해독·부가 프로세스를 제외한 약1/3을 sort에 배분함. 68CPU·두HAP이면 각각 정렬17/sort8/입력해독4/부가5로 합계68임. 숫자 명시값·RAM 모델·sort 버퍼 축소 규칙은 유지함. 이미 실행 중인 도구의 스레드는 변경하지 않음.
- 독립 source trace : writer와 별도로 source 필터·CIGAR 투영·중복 처리·전체 SHA256 multiset·depth를 검증함. canonical JSON 바이트 규약을 보존하는 C++/OpenSSL 배치 해시를 사용함. source는 최대4096 records 또는16MiB에서 배치를 비우며, Zarr는 기존 packed16MiB 상한을 유지함. read와 segment의 JSON/SHA256을 배치 처리하고 multiset을 연속 기록함. contig/HAP 병렬 worker 안에서 CPU 하나로 수행함. SHA context를 배치 안에서 재사용함. `MFM_TRACE_HASH_BACKEND=auto|native|python`이며 auto 빌드 불가 시 Python을 사용함. 실제 backend를 QC에 기록함.
- 환경 : HTSlib·OpenSSL·C++17 컴파일러를 conda `environment.yml`과 00 설치의 버전 기록 대상에 포함함. pip 전용 의존성은 늘리지 않음. 실행 중 환경은 재설치하지 않으며 현재 jj_mfm에서 실제 빌드를 검증함.
- 입력 재사용 : `MFM_GPU_INPUT_CACHE_ROOT`를 지정하면 같은 완료 지문의 FASTQ 캐시를 별도 비교 출력에서도 공유 가능함. 키·paired 여부·파일 크기/mtime를 기존과 같이 검증하며 입력이 다르면 재사용하지 않음.
- 검증 범위 : 실제 작은 BAM에서 직렬/병렬 CRAM의 모든 SAM 필드·순서·41 records가 일치함. native SHA는 기존 JSON SHA와 일치하고 입력 변조를 검출함. 실제 NA18959 chr22 HAP1 전체121,614 retained alignments의 필터·multiset·segment·depth 결과가 일치함. 전체 trace wall은61.79→61.75초로 유의한 시간 단축은 아직 확인되지 않았으며 CPU 시간은36.29→34.23초임. 추가한 read/segment 배치 최종본도 같은 실제 구간의 전수 결과 일치를 확인함. 세부 wall은 공유 부하의 영향을 받으며 전체 trace 가속률로 단정하지 않음. 두HAP 중첩의 전체 샘플 가속률은 별도 비교 중임.

- 전량 CRAM 실측 : HG00275 HAP1의 동일675,604,408 records BAM에서 직렬389.72초→13worker 병렬265.70초로 관측1.47배임. 전체 비교 대상 read 필드 checksum 및66개×1024bp depth 구간이 일치함. 전체 genome 모든1bp를 별도로 대조한 결과나 CRAM byte 동일 검사와 구분함. 검증 후 이 시험이 만든 BAM 링크·CRAM·shard를 정리함.
- 4GPU Illumina 옵션 실측 : 같은 전량 FASTQ/개인 HAP index에서 기본464.78초, CIGAR만501.12초, partition322.79초, 두 옵션251.31초였음. 네 조건 모두 전체 read checksum·depth probe를 통과함. 공유 CPU 부하가 달라 보편적인 최적값이나 전체 전처리 가속률을 확정한 결과는 아님. 4GPU 검증 profile로 사용하며 명시 설정을 우선함.

## RNA 의존성 병렬화·네이티브 I/O (2026-10-31)

- 실행: 기존 `run_pipeline.sh run/resume`에서 내장 `run_rna_sample.sh`가 수행함. 별도 운영 제어기 불필요함. 이미 실행 중인 계산에는 소급 적용하지 않음.
- HAP 분기: 각 HAP가 `15 index → 16 align → 17 duplicate 표시 → 18 coverage → 19 HAP QC`로 진행함. HAP1 index가 먼저 끝나면 HAP2 index 종료를 기다리지 않고 HAP1 정렬에 진입함.
- WGS 기반 분기: `20 variant 검증`과 `22 CNV`는 두 HAP STAR 과정과 겹쳐 진행함. 두 `17`과 `20`이 모두 성공하면 `21 ASE`를 coverage와 겹쳐 진행함. HAP QC·variant·ASE·CNV가 모두 통과한 뒤 summary·23 ASE QC·24 metadata·25 최종 QC·26 cleanup 진행함.
- 재개/중단: 명시한 `RNA_START_AT>=15` 또는 `RNA_STOP_AFTER`, index benchmark, HAP jobs=1, 작은 CPU 예산은 번호별 barrier 경로를 사용함. 실패하면 후속 최종 검증을 차단하고 유효한 형제 분기 결과는 보존함. 종료 신호는 해당 분기 자식·CPU/GPU lease까지 정리함.
- CPU: affinity/cgroup/Slurm의 작은 한도를 사용함. 현재 cgroup68코어임. 풀을 줄이지 않고 동시 분기별 상한을 나누며, 고정 도구의 스레드는 종료 경계에서 재산정함. 예약 예외 GPU 비교도 분기 상한을 따름. `qc/rna_dag_execution.json`은 CPU 상한과 시작/종료 시각이며 실제 사용량은 단계 resource audit의 CPU초/경과시간으로 구분함.

| 설정 | 현재 기본/동작 |
|---|---|
| `RNA_PIPELINE_OVERLAP` | `auto`; 위 선행조건에서 의존성 병렬 실행함. `0`은 번호별 barrier임. `1`도 명시 STOP/재개·자원 조건을 우선함. |
| `RNA_STAR_INDEX_HAP_JOBS` | `2`; 명시한 `1`은 보존함. |
| `RNA_ALIGN_HAP_JOBS` / `RNA_COVERAGE_HAP_JOBS` | 각각 `2`; 번호별 경로에서도 같은 총 CPU 예산을 나눔. |
| `RNA_MARKDUP_HAP_JOBS` | 기존 `2` 유지함. |
| `RNA_ALLELIC_JOBS` | `auto`; 실제 단계 CPU와 RAM 예측 안에서 정함. 명시 양수는 상한임. |
| `RNA_ALLELIC_INPUT_BACKEND` | `auto`; 합산 BAM64MiB 이상이면 네이티브 QNAME shard 사용함. `partitioned`/`legacy`로 명시 비교 가능함. |
| `RNA_DEPTH_BACKEND` | `auto`; HTSlib C++로 contig를 한 번 읽음. `native`는 빌드 실패도 오류임. `legacy`는 기존 samtools/구GPU 경로임. |

- Coverage: HTSlib가 MAPQ255·primary/QC pass·duplicate 포함·CIGAR D 포함/N 제외·mate overlap2회 계수를 그대로 적용함. CPU는 binary uint32, GPU는 같은 int64 차분의 누적합을 계산함. GPU 자식은 계산·전달 후 종료하고 CPU Zstd5 쓰기 동안 CUDA context/장치 예약을 보유하지 않음. 작은 window QC는 CPU에서 계산함.
- 전장/window 저장: 같은 임시 raw buffer에서 두 Zarr를 작성하여 전장 Zarr의 반복 읽기·해제를 줄임. flat chunk당 한 번만 압축하고 반환한 CPU를 실제 단계 상한 안에서 압축 전에 재확보함. 기존 배열·행 순서·dtype·raw depth SHA를 유지함. buffer 크기·mtime·SHA를 검사하고 종료 시 자동 제거함. 최대 contig 길이·worker 수로 RAM과 임시공간을 계산하며 작은 실제 CPU/RAM 할당을 우선함.
- ASE 입력: 양 BAM을 QNAME hash의 같은 binary shard에 한 번씩 기록함. BGZF 해독 thread·입력 job·부모 CPU를 단계 grant 안에 배분함. 모든 mate/교차염색체/양 HAP 증거가 같은 shard에 있으므로 duplicate 합집합·AS 동점·길이 보정·공통 gene·검증 variant·split-half 규칙을 보존함. 작은 count만 합치고 임시 shard는 종료 시 제거함. 입장 CPU가 작아도 큰 입력을 legacy로 바꾸지 않으며 단계 상한 안에서 CPU를 재확보함. fork worker의 RAM은 join까지 유지하고 tmpfs 임시파일 예상량도 예약에 포함함. worker 상한과 실제 overlap/CPU초를 구분해 기록함. beta-binomial/BH FDR는 기존23에서 계산함.
- FASTQ/중복 표시: gzip 해제·mate suffix 복원·count/최대 길이를 한 번에 계산하며 mate별 동시 처리함. Parabricks `bam2fq` 자체는 CPU 도구임. 분리 경로의 필수 queryname sort와 통합/분리 모두의 duplicate/optical 기준·metrics를 보존하며 독립 flag QC의 Python loop를 HTSlib 단회 집계로 대체함.
- 환경: 기존 `environment.yml`의 HTSlib·C++ compiler와 zlib 의존성 사용함. 네이티브 도구는 기존 빌드 캐시에 저장함. auto의 컴파일 불가에만 기존 경로로 복귀하며 손상된 입력을 fallback으로 숨기지 않음. STAR GTEx filtering과 GPU strict/experimental 정책은 변경하지 않음. 대형 GPU 차분은 예측 VRAM mapping만큼 CUDA 자식의 가상 주소 공간을 추가하며 실제 CPU RSS 감시는 유지함.
- 검증: D/N·MAPQ·duplicate·mate overlap·빈 contig·FASTQ 손상·QNAME 동점·모든 ASE TSV·분기 실패/신호/명시 STOP를 대조함. 완료 HG00126의 동일 RNA BAM15.70GB로 전량 ASE와 HAP1 전장/window depth를 별도 경로에서 재계산함. 아래 실측은 해당 단계의 가속률이며 샘플 전체 가속률로 확대하지 않음.

동일 HG00126 전량 실측(2026-10-31, 큐 대기 제외):

| 단계/경로 | 기존 | 개선 | 단계 속도 향상 |
|---|---:|---:|---:|
| ASE 양 HAP BAM15.70GB | 689.878초 | 49.054초 | 14.06배 |
| HAP1 전장+2,640창 coverage/CPU | 74.296초 | 25.595초 | 2.90배 |
| 같은 coverage/GPU | — | 52.766초 | 개선 CPU보다 느림 |

- ASE4개 TSV는 기존 완료본·legacy 재실행·최신 경로 모두 SHA256 일치함. depth는 chr1–22 및 window store의 모든 배열·값·shape를 기존 완료본/legacy/최신 CPU/GPU에서 대조해 일치함.
- CPU 상한은 각68이지만 평균 실사용은 ASE1.52→20.64코어, coverage11.06→3.02코어임. coverage CPU 누계822.07→77.29초로 계산량 자체가 줄었으므로 낮은 평균 코어가 곧 성능 저하를 뜻하지 않음.
- 공통 chunk writer 개선을 양쪽 기준선에 적용한 단계 비교임. 동시에 다른 작업이 있던 단일 샘플 관측이며 전체 RNA/다운로드/STAR 총 가속률과 구분함. GPU context는 prefix 자식 종료 후 반환되고 CPU 압축 중GPU0–3 메모리0을 확인함. 기본 CPU 선택은 유지하며 명시 CuPy 선택도 지원함.
- 원본/전달용 각각42개 관련 검사 통과함(GPU 선택 검사는 별도 실제 실행함). 실제 Parabricks/Picard duplicate flag·read·필수 태그 대조도 통과함. PG는 도구 실행 이력 태그이므로 비교에서 분리하며 나머지 태그와 read 필드 전체를 대조함.

## HAP 변이 검증·CNV의 WGS 읽기 가속 (2026-10-31)

- 적용: 원본20/22, 전달본18/20에 연결함. 기존 `run_pipeline.sh`의 RNA 실행 흐름을 사용하며 별도 운영 제어기는 추가하지 않음. 기존 실행과 완료 결과는 보존함.
- 변이 지원도: VCF 영역 선택·chain 변환·HAP allele 확인·gene 부여·판정 코드는 유지함. 각 염색체의 조회 범위를 HTSlib로 한 번 순회하고 CIGAR block에서 SNV/MNV·indel 지지수를 계산함. 긴DEL이 여러 원래 영역에 조회되는 occurrence도 각각 유지함.
- 변이 필터: MAPQ20 기본·unmapped/secondary/supplementary/QC-fail/duplicate 제외함. HP별지지수·indel종류/길이·위치±10bp·양쪽flank·splice N제외·최소깊이10/alt0.1 판정 유지함. 기존pysam의P offset호환·잘못된float HP1/2실패도 보존함.
- CNV: exon 경계점에서 read span의 겹침 길이를 누적함. 전장1bp배열을 만들지 않고 BED각행의합을 반환함. MAPQ0·UNMAP/SECONDARY/QCFAIL/DUP 제외·supplementary포함·D/N포함·mate겹침2회 계수임. 설치된samtools1.24와 직접 대조함. [공식 bedcov 설명](https://www.htslib.org/doc/samtools-bedcov.html)의 기본 규칙도 확인함.
- CNV 정책: gene별exon합집합·HAP별전체gene평균의중앙값·최소길이100bp·결측처리·RD0.5–1.5·두HAP판정·row순서·TSV/metadata schema 유지함. RNA coverage의MAPQ255/duplicate포함 필터와 별개임.
- 자원: 실제CPU68·RAM1.65TB/예약1.485TB와더작은실제할당 안에서 실행함. 68CPU에서 최대22작업×(계산1+HTSlib해독1+I/O1)+부모1=67CPU 계획임. worker4GiB·부모4GiB의상한예측을먼저기록함. 변이fork worker RAM은join까지 유지함. CPU는작업경계에서반환/재확보하고DAG분기상한은 유지함. CPU할당과실사용은audit에서구분함.

| 선택값 | 동작 |
|---|---|
| `RNA_HAP_VARIANT_BACKEND=auto` | 지원도C++ 도구를 사용함. `native`는빌드실패도오류, `legacy`는기존Python경로임. |
| `RNA_CNV_BACKEND=auto` | 희소구간합C++ 도구를 사용함. `native`/`legacy(samtools bedcov)`로명시비교 가능함. |

- 환경/실패: 기존HTSlib·C++컴파일러 및 native cache를 사용함. auto는컴파일불가일때만legacy임. 잘린입력·index/reference누락·잘못된값은실패하며숨기지않음. 입력/출력 임시TSV는종료시정리함. Native CNV는contig당retained INT32_MAX초과 시bedcov depth cap의정확성검증을요구하며실패함.
- 검증: SNV/MNV·indel·softclip/P/I/D/N/=/X·HP·영역경계/중복occurrence·BAM/CRAM·exon겹침·유전자명공백·빈chr·index누락·도구실패를legacy/samtools와대조함.

HG00126 전량 비교 결과:

| 단계 | 과거 완료 로그 | 이번 native 실행(큐 제외) | 출력 대조 |
|---|---:|---:|---|
| HAP 변이741,235행 | 29분11초 | 69.539초 | 전체TSV SHA·QC 일치함 |
| CNV75,035 gene | 18분9초 | 64.869초 | 전체TSV SHA·정책/중앙값 일치함 |

- 위두실행은서버부하/시점이다르므로순수알고리즘가속률로해석하지않음. 전체RNA13–26의새벽시계시간은이번측정에포함하지않음.
- 추가동조건검증: 같은chr22·CPU상한4에서변이검증legacy274.677초/native21.204초, 평균실사용각0.997/1.004코어였으며TSV전체SHA가같음. 이chr22직접비교에서는약12.95배임. 전량가속률로확대하지않음.

- 실제chr22/HAP1 CRAM·같은BED9,298구간·CPU상한4 직접비교에서도bedcov 150.149초/native 5.950초였고원시구간합TSV SHA가일치함. 이구간비교는약25.24배이며전량가속률과구분함.

## GPU RNA 정렬·중복 표시 통합 — 2026-10-31

- 내장 `run_pipeline.sh`/`run_rna_sample.sh`의 `RNA_GPU_MARKDUP_MODE` 기본값은 `integrated`임. RNA와 BAM backend가 모두 Parabricks일 때 적용함. CPU STAR·Picard 경로의 정렬/필터는 유지함.
- 16번: 기존 STAR filtering을 유지한 `rna_fq2bam` 안에서 좌표정렬·중복표시·metrics를 함께 생성함. 설치된4.8.0-1의 optical 거리100 인자와 결과를 실제 검증함. queryname-assume 옵션은 사용하지 않음. GPU STAR의 기존 `Normal`/CPU GTEx의 `BySJout` 차이는 그대로 명시하며 strict 정책은 계속 거부함.
- 17번: 통합 BAM/metrics/adapter fingerprint·기존 검증표식을 확인한 뒤 독립 BAM 복사·인덱스·네이티브 flag QC를 수행함. 원본 align inode에 hardlink를 추가하지 않아 cache의ctime을 보존함. 최종 출력과 나머지 RNA 단계는 기존 계약을 사용함.
- 모드·adapter·metrics 변경은 정렬 cache를 무효화함. 비교/복구 시 JSON `pipeline_env`에 `"RNA_GPU_MARKDUP_MODE": "separate"`를 지정 가능함. 실행 중 변경은 다음 정렬부터 적용하며 유효한 기존 계산을 종료하지 않음. 직접 STAR adapter 호출은 기존 분리 동작을 유지함.
- 자원: GPU RNA runner의 host RAM 기본값은HAP당256GiB이며명시한 `MFM_GPU_HOST_MEMORY_BYTES`를우선함. 실제 CPU quota/Slurm/단계 grant를 준수함. 직접 MarkDuplicates Python 실행도cgroup quota를반영함. 실물시험은CPU상한68·두HAP각34thread·GPU0–3중각1장·host상한각256GiB를사용함. 전체RAM1.65TB/예약1.485TB와별도프로젝트불변원칙을유지함.
- 실물 검증: HG00126 두HAP의220,192,091개정렬레코드(FLAG·서열/품질·이름·좌표·MAPQ/CIGAR·mate·PG외모든tag), duplicate/optical metrics, 전장/window Zarr의모든배열, ASE·CNV등7TSV를기존GPU완료본과대조PASS임. 작은광학경계시험88096레코드는Picard와전필드/metrics동등임. 최종105QC PASS임.
- 관측시간: 두HAP 정렬+중복표시의성공실행구간 합집합은과거46분 51초→이번9분 10초임. 누적개선을포함한RNA16–24는과거136분 27초→이번12분 30초임. 같은입력이지만부하·thread배분이달라관측값으로보고함.
- 측정 범위: FASTQ/lift/STAR index와WGS를재사용함. 이번속도시험은CPU native coverage를명시선택했으며GPU mode의coverage backend기본선택은바꾸지않음. 최초준비·재검증·시험사본복제시간은계산구간과분리함. 26정리는benchmark보존을위해0으로두었으며일반성공후cleanup정책은유지함. 새cold 전체RNA가속률은미측정임.

## RNA 후속 CPU·FASTQ 보완 — 2026-11-01

- 동적 CPU: 내장 DAG의 변이(20)·CNV(22)·ASE(21)는 처음1/4풀로 입장하며 maximum/throughput·본인RNA공유lease·정상예약이 확인된 native 작업만 작업 경계에서 CPU를 재확보함. 구간 모드는 부모1+작업당1CPU, contig 대조 모드는 부모1+작업당3CPU(계산/해독/I/O)임. 다른 lease·RAM·Slurm·명시 상한을 지키며 실제 사용코어와 예약 상한을 구분함.
- 고정 STAR/Picard·명시jobs·balanced·비예약비교·대여상한은유지함. CPU/RAM은실제in-flight worker보다먼저반환하지않으며fork RAM은join까지보존함. DAG JSON의 `cpu_cap`은입장상한, `elastic_cpu_upper`는조건부재확보상한임. 실제grant/peak/CPU초는CPU원장으로확인함.
- FASTQ(13): `RNA_FASTQ_BACKEND=auto|samtools|parabricks`이며 `auto`는 검증된 samtools 비압축 경로를 사용함. GPU 정렬·중복표시는 유지함. JSON `pipeline_env`에서 `"RNA_FASTQ_BACKEND": "parabricks"`를 명시해 PB 경로를 선택 가능함. 요청·실제 backend·cache/provenance를 기록함.
- samtools pipe는실제grant 안에서두계산·BGZF해독/I/O thread와예상RAM을먼저산정함. PB/samtools모두변환종료후네이티브검사단계에서CPU를필요한worker수+부모1로낮춤(5개출력검사이면최대6). 검사완료후1CPU와부모RAM만남김. 고정변환도구실행중thread수는바꾸지않음.
- 검증: 같은 원본의 HG00126 FASTQ 36,602,796pair를 mate별 이름·서열·품질 checksum으로 대조함. 출력 순서 차이는 있으나 두HAP GPU 후속 BAM의 PG 외 모든 필드/tag·metrics, 4Zarr 전배열, 7TSV까지 전량 대조 PASS임.
- 기존 PAIR2 보호 예약을 유지하고 비교 종료 후 공식 공유풀에서 RNA 후속 시험을 완료함. GPU0–3 범위·RAM1.65TB/예약1.485TB를 지켰으며 GLM·기존 계산은 변경하지 않음. WGS/lift/index는 재사용했으므로 새 cold 전체 전처리 가속률로 해석하지 않음.

### 동조건 전량 검증 후 적용한 추가 개선

- **FASTQ 준비(13)**: 동일68CPU·64GiB 예약의 PB→SAM→SAM→PB 비교에서 PB223.9/224.2초, SAM97.2/98.9초였음. SAM 기본8/4 thread에서도105.3초·read집합 동등을 확인함. config의 기본 collate/fastq는8/4이며 사용자 명시값·실제 grant로 제한함. auto 기본값을 SAM으로 변경함.
- **예약 반환**: 정상 RNA 변환 lease에서 samtools 실행 전 실제 고정 thread 상한·예상 RAM으로 예약을 낮춤. 기본8/4는 부모/해독/I/O 포함 최대16CPU이며 남는 CPU/RAM을 다른 단계에 반환함. 실행 중 thread 수는 바꾸지 않으며 대여·비예약·BENCH 조건은 유지함.
- **변이(20)·CNV(22)**: `RNA_WGS_TASK_LAYOUT=auto|contig|chunks`, `RNA_WGS_CHUNK_BP=16000000`을 사용함. auto는 CPU풀8 이상에서 chunks를 선택함. 원래 병합 영역/BED행은 자르지 않고 원래 occurrence·행 번호·합산 순서를 보존함. native 추가 해독0·작업당1CPU·부모1CPU로 최대67worker/68CPU 안에서 실행함. 다른 lease의 RAM을 제외해 worker 상한을 계산하며 worker당4GiB·부모4GiB 예산 및 주소공간 제한을 유지함. 명시 threads/jobs·단계 상한은 유지함.
- **구간 실측**: 같은68CPU 예약에서 변이 contig55.8–70.1초/16Mb29.3초, CNV contig65.8–65.9초/16Mb33.2초였음. 8Mb는28.8/35.5초로 더 작게 나눌 일관된 이득이 없어16Mb를 기본으로 선택함. 모든 전량 TSV SHA가 같음. task별 span/행수/벽시계/CPU초를 QC metadata에 기록함.
- **ASE(21)**: `RNA_ASE_SHARD_READER=auto|python|native`이며 auto는 C++/Python C API reader를 사용함. 기존 jj_mfm의 C++컴파일러·Python.h로 빌드하며 새 pip 의존성은 추가하지 않음. 빌드 불가 시에만 Python으로 전환하고 손상 데이터 오류는 숨기지 않음. Unicode·AS동점·mate·duplicate·MAPQ·입력 순서 및 모든 count/meta를 전량 검증함. 해석/필터/배정 CPU시간을 분리 기록함.
- **ASE 실측**: Python40.997/41.490초·CPU1000.6/1023.7초, C++40.453/38.842초·CPU909.7/912.8초였음. 시간 이득은 작으며 약10% CPU 누계 절감이 주효과임. 모든 count TSV와 생물학적 metadata가 같음.
- **해석 범위**: 동일 샘플의 단계별 관측값임. OS page cache·다른 cgroup 작업 부하는 통제하지 않았고 단계 간 overlap이 있어 절감 시간을 단순 합산하지 않음. 공식 내장 .sh 실행에 적용되며 상시 감시기·별도 샘플 운영큐는 추가하지 않음.


### 2026-11-01 GPU 입력 공급·Short 준비

- WGS FASTQ 검사·mate 분할은 `MFM_GPU_FASTQ_BACKEND=auto`에서 paired 입력에서 CPU와 빌드 환경이 맞으면 네이티브 경로를 사용한다. `python`은 기존 대조 경로, `native`는 명시 실행이다. 행 상한64MiB·프레임4MiB·형식/QNAME/길이 검증과 FASTQ 바이트·태그를 보존한다. 공통 캐시 완료 장부에 backend와 압축 CPU 배분을 기록한다.
- Short의 `MFM_GPU_INPUT_MODE=auto|cache`는 `MFM_GPU_PREPARE_OVERLAP=1` 기본에서 같은 BWA index 생성과 FASTQ 캐시 준비를 겹친다. index1CPU와 입력 CPU/RAM/디스크를 합산 예약하며 캐시 잠금 대기에는 index 자원만 유지한다. CPU/RAM이 작거나 저장 파일시스템이 다르면 순차 경로를 쓴다. `stream` 명시 또는 overlap0은 기존 준비 순서를 유지한다. 인덱스 알고리즘/내용은 변경하지 않는다.
- 준비 중첩 시간·세부 CPU/RAM은 GPU 정렬 metadata의 `input_preparation_overlap`에 기록한다. 이때 `input_seconds`에는 index와 입력의 중복 제거 경과가 포함되므로 기존 입력 단독 시간과 바로 비교하지 않는다.
- RNA는 기존 `RNA_FASTQ_BACKEND=samtools`와 `MFM_COVERAGE_BACKEND=cpu`를 GPU 중복 표시와 조합할 수 있다. FASTQ 변환 도구와 STAR/중복표시/coverage backend를 별도로 선택한다. 전체 속도는 캐시 최초 생성과 후속 결과 동등성까지 측정해 판단한다.

- Long은 부분집합 AB에서 네이티브 경로의 이득이 없어 auto가 기존 Python 경로를 사용한다. FIFO 공급 CPU는 `MFM_GPU_STREAM_CPUS`로 명시할 수 있으며 기본16 이내다. 실행 중 공급 CPU/RAM은 GPU 대기 중에도 반환하지 않고 계산 grant에 포함한다. 32CPU 확대는 전량 시험으로 확인 후 결정한다.

- **자원 절충**: CNV16Mb는 CPU누계845→1167초(약38% 증가), RSS합 표본최대21.8→62.3GB를 사용해 벽시계66→33초로 줄었음. 변이의 RSS합은14.4→39.5GB였음. 작업당 주소공간 제한·공유RAM 예약을 유지함. 전체 동시샘플 처리량 최적성은 별도 미측정임.


### 2026-11-02 RNA 통합 가속·임시 RAM 예약 보완

- 적용: ASE 입력 분할 후 실제 spool의 논리/할당 바이트 중 큰 값을 RAM 예약에 반영함. `/tmp`는 tmpfs이므로 RAM 집계에서 제외하지 않음. 모든 worker 종료와 임시 폴더 정리 후 부모 RSS+여유만 예약하고 CPU를 반환함. 공유 풀의 다른 작업 예약을 포함해 worker 수를 계산함.
- CNV: `RNA_CNV_CRAM_FIELDS=full|coverage` 지원함. `coverage`는 FLAG/좌표/CIGAR/SEQ를 해독하고 QUAL/aux 등 미사용 필드를 생략함. SEQ·참조 MD5 검증을 보존하며 잘못된 reference는 실패함. 전량 동등성·성능 검증을 통과하여 기본값을 `coverage`로 확정함.
- ASE: `RNA_ASE_ASSIGN_BACKEND=tables|auto|native` 지원함. C++에서 shard 읽기→중복 처리→길이 보정→HAP/유전자 배정을 연결해 Python read 객체 생성을 줄임. split-half는 기존 Python 해시 함수를 호출함. 기존 21/19번의 배정 규칙과 함수 지문이 일치할 때만 사용함. `auto`는 미지원 정책/범위에서 기준 경로로 전환하며, 손상 입력은 실패함. 전량 동등성·성능 검증을 통과하여 기본값을 `auto`로 확정함. 검증된 정책에서 native를 우선 사용하며 미지원 정책/빌드 불가 시 tables를 유지함.
- 유지: MAPQ·duplicate·mate·AS 동점·evidence·window/N/blacklist·RNA 신호 및 ASE 생물학적 기준을 변경하지 않음. 실행 옵션은 기존 JSON의 `pipeline_env`에도 지정 가능함. 추가 운영 제어 스크립트 없이 `run_pipeline.sh run --selection ... --config ... --executor local` 경로에 적용됨.
- 검증: 무작위/경계·잘못된 참조·실제 BAM CLI·임시 예약 반환 검사 통과함. 전량 ABBA 대조 후 TSV/metadata 일치와 시간5% 또는 CPU10% 개선(시간 회귀5% 이내)이 확인된 후보만 기본값으로 승격함. 전체 속도 향상 수치는 별도 단일/2샘플 병행 시험으로 판단함.


### 2026-11-02 실측 유휴 CPU의 논리 대여

- `cpu_backfill.max_claim_pools`는명시정책에서최대4, `max_loan_pools`는최대3까지지원한다. quota68환경의272논리예약/204대여한도는실제272코어를할당받는뜻이아니다. 기본값2/1과기존3/2설정은그대로다.
- 두입장경로는같은`loan_limits`계산을사용한다. 10초실제유휴관측·대여후안정화·실제cgroup quota일치·RAM예약/cgroup/RSS보호를통과해야새고정단계를대여한다. 실행중worker의thread/nice/affinity와입력/QC는변경하지않는다.

- 명시 `cpu_backfill.preferred_samples` 목록이 있으면 해당 샘플의 대여 입장을 먼저 검토하고 같은 우선순위 안에서는 기존 FIFO를 유지한다. 목록이 없으면 기존 순서다.

최종 BAM 쓰기·전수 집계가 끝난 색인/QC 단계는 parent 4GiB와 추가 스레드를 포함한 스레드당 64MiB로 RAM 예약을 반환합니다. 이전 쓰기/읽기 검증의 예약 및 CPU 수는 유지합니다.


### 2026-11-02 GPU BWA 인덱스 재사용·선행 준비

- 인덱스 식별: FASTA 전체 바이트 SHA256·BWA 실행파일 SHA256·생성 옵션을 사용함. 파일 경로/mtime이 달라도 같은 입력이면 같은 캐시 경로에서 재사용 가능함. 서로 다른 개인 HAP나 다른 BWA 실행파일은 공유하지 않음. BWA-MEM2 인덱스와 혼용하지 않음.
- 무결성: FASTA·FAI·BWA 구성파일5개의 크기/SHA256 및 완료표식을 검사함. 잠금으로 중복 생성을 막고 임시 폴더에서 생성·검사한 뒤 원자적으로 공개함. 생성 실패·입력 변경·부분 결과·손상 캐시는 사용하지 않으며 기존 캐시를 자동 삭제/덮어쓰지 않음.
- 준비 중첩: `MFM_GPU_PREPARE_OVERLAP=1`이 기본임. GPU 입력 BAM이 준비된 뒤 인덱스가 없으면 인덱스1CPU와 FASTQ 최초 생성을 합산 CPU/RAM/디스크 예약 안에서 겹침. `stream` 요청도 최초 인덱스 생성 시에는 FASTQ 임시 캐시를 만들어 중첩함. 인덱스가 있으면 원래 stream/cache 경로를 유지함. CPU/RAM이 부족하거나 저장장치가 다르면 순차 경로를 유지함. upstream BAM 준비와의 중첩은 이번 적용 범위가 아님.
- 자원 반환: 입력 생산·압축 프로세스가 모두 종료하면 입력 CPU 몫을 반환하고 인덱스1CPU·예상 RAM만 유지함. GPU를 확보할 때 계산 자원을 다시 예약함. 최초 인덱스 계산 자체는 CPU 직렬 구간임.
- 보존 범위: 기본 캐시는 해당 GPU work의 `bwa_indexes/`임. 재시도·동일 작업 경로에서 재사용함. 실행 간 재사용이 필요하면 기존 JSON `pipeline_env`에 `MFM_GPU_BWA_INDEX_ROOT`로 공통 캐시 경로를 지정함. 명시한 외부 공유 캐시는 모델 패키지 정리로 삭제하지 않으므로 필요한 비교/재시도 기간만 보존하도록 운영함.
- 기록: HAP별 GPU JSON의 `bwa_index_cache`에 생성/재사용 여부·서열/도구 지문·준비시간·독립 인덱스 CPU/RSS 측정을 저장함. `input_preparation_overlap`에 인덱스/입력 준비시간 및 입력 완료 뒤 예약1CPU를 기록함. 정렬 명령·read/tag/duplicate·필터 및 최종 모델 입력 정책은 동일함.


### 2026-11-02 개인 인덱스와 upstream BAM 준비 중첩

- 기본 활성: `MFM_GPU_BAM_INDEX_OVERLAP=1`임. 기존 `run_pipeline.sh`/`run_sample.sh` 안에서 내부 worker를 실행하므로 별도 운영 제어기가 필요하지 않음. JSON `pipeline_env`에 `0`을 지정하면 기존 순차 준비를 재현함.
- 적용: 단일 샘플의 GPU Illumina·all_to_both·입력검증 활성 경로임. 개인 target FASTA를 먼저 확정하고, 아직 없는 HAP 인덱스와 routed/name-collated BAM 준비를 병행함. HAP 인덱스가 모두 준비돼 있으면 기존 BAM 준비 경로를 사용함. HiFi/다른 routing/dry-run은 기존 경로를 유지함.
- 예약: HAP별 인덱스1CPU+부모제어1CPU를 기초몫으로 두고, 남은 CPU·CPU당2GiB의 입력 예산과 인덱스 RAM 추정치를 같은 lease에 합산함. 예를 들어 전체 grant68이면 입력65+인덱스2+제어1 이내임. 입력 worker의 routing/해독/집계/압축과 주소 공간은 해당 하위 예산으로 제한함. 전체 CPU풀을 줄이지 않으며 실제 grant·Slurm/cgroup 상한을 우선함.
- 반환: 자식과 하위 프로세스의 종료를 확인한 뒤 해당 몫만 반환함. 입력이 먼저 종료되면 인덱스2+제어1을 유지하며, 인덱스가 먼저 종료되면 입력+제어만 유지함. 두 준비 완료 후 GPU 계산 자원을 다시 확보함. 실패 시 이 호출의 자식 그룹만 회수함.
- 저장공간: 동시 BAM·인덱스 쓰기 예상량을 장부에 기록하며 다른 GPU 작업의 준비 예약도 합산함. RAM/디스크가 부족하거나 서로 다른 저장장치이면 순차 경로를 사용함. 기존 캐시·원본을 이동/삭제하지 않음.
- 정확성·기록: read 선택·태그·duplicate·검증 규칙과 BWA 생성 명령을 유지함. `logs/personal_realign/<sample>.prepare_overlap.json` 및 stage timing에 각 worker의 실제 CPU초/시간·RSS 표본·중첩시간 합집합을 기록함. RSS 합계는 공유 페이지를 중복 계상할 수 있는 프로세스 트리 표본이며 CPU 예약과 실제 사용량을 구분함.

- 선행 다운로드: 명시한 `download_allowed_block_reasons`에 한해 계산 보류와 다운로드 보류를 분리합니다. 일반 준비 버퍼가 차도 미준비 pair 원본은 우선 확보하며, `prefetch_max_bytes`로 전송 중 입력 통계 합계와 디스크 안전 여유를 확인합니다. 큐의 전송량 변경은 실행 중 샘플의 데이터·계산 계약을 바꾸지 않습니다.

- 전송 RAM은 기본1GB이며, Python/스레드512MiB와 HTTP 범위worker별3개32MiB버퍼에20%여유를합산해병렬범위수가크면자동증가합니다. checksum과범위/원격동일성검증은유지합니다.

### 2026-11-02 검증 작업의 RAM 입장·대기 반환

- 예측 대기: `WorkerBudget(..., wait_memory_bytes=...)`를 명시한 계산 전 대기는 미할당 예측량을 유지하지 않음. 실제 부모 RSS+10%와 지정한 보존 버퍼 최소량을 남김. 실제 worker가 실행 중이면 CPU/RAM을 반환하지 않음. 기존 호출의 기본값은 `None`으로 유지하며 fork worker·spool 등 보존 자료가 있는 호출에 일괄 적용하지 않음.
- 검증 시작: `BENCH`만 `MFM_BENCH_ENTRY_MEMORY_BYTES`로 시작 RAM을 명시할 수 있음(최소4GiB, 미지정64GiB 유지). 실제 RNA/WGS 단계의 시작 예약에는 영향을 주지 않음. 큰 계산을4GiB로 수행한다는 의미가 아니며 계산 전 필요한 RAM을 다시 확보해야 함.
- 전량 비교 운영: CNV는 context·부모·native worker별4GiB, ASE는 입력 BAM 크기·분할 spool·해독 버퍼·assignment worker를 반영해 RAM을 예측함. 가용 공유 RAM 안에서 동시 worker를 정하며 각 ABBA4회에는 같은 단계 한도를 유지함. 단계가 모두 종료되고 임시 자료가 정리된 뒤4GiB로 반환함. 대기·예약·실제CPU/RSS는 별도로 기록함. CPU 상한68과 운영RAM1.65TB/예약1.485TB를 유지함.
- 적용 범위: 운영 전처리는 기존 `run_pipeline.sh`를 그대로 사용함. 별도 전량 검증의 자원 대기를 개선했으며 read/window/태그·출력 규칙은 바꾸지 않음. 전체 가속률은 전량 비교 완료 뒤 판단함.

- 전량 검증 후 기본값 반영(2026-11-02T13:11): CNV CRAM필드=coverage, ASE선택=auto(검증된 정책에서 native 우선)임. gene_rd/ASE4TSV·생물학적metadata 동등 및 승격 후 회귀검증을 통과함. 전체·병행 시험은 별도 진행함.

- 전량 확정 비교(2026-11-02): HG00126의 같은 입력으로 ABBA 각4회/방식별2회 평균을 비교함. CNV는 같은29GiB·6worker에서172.45→55.29초(3.12배), ASE는 같은113GiB·22assignment worker에서85.84→61.13초(1.40배)임. CPU 누계는 각각70.27%·72.63% 감소함. CNV 전체TSV와 ASE4TSV의20개 대조 및 ASE 생물학 metadata가 모두 같음. 실행 CPU상한68·16Mb를 유지하며 RAM대기시간을 제외함. 외부부하/page cache 미통제·단일샘플 결과이며 전체 파이프라인 가속률은 별도 확인함.


### Short WGS–RNA pair의 DNA 서열만 먼저 처리함

- JSON `pipeline_env`에서 `MFM_SHORT_READ_OUTPUT_MODE=sequence_only`, `MFM_SEQUENCE_DEFER_RNA=1`을 함께 지정함. 기본값은 RNA 보류하지 않음(`0`)이며 Long에는 미적용임.
- `run_pipeline.sh run --selection <기존 pair 선택 TSV> --config <설정 JSON> --executor local`로 실행함. 기존 `data_group=wgs_rna_pair`와 `rna_path`/URL을 보존하며 별도 output/state 경로를 사용함.
- 공통 phased VCF/reference로 HAP1/HAP2 서열을 생성하고 기존 window 정책·서열 Zarr·좌표매핑·독립 source trace·최종 검증을 적용함. WGS/RNA read 다운로드·정렬·coverage·ASE·CNV는 실행하지 않음.
- 최종 검증의 `mode=wgs`, `completion_scope=dna_sequence_only_rna_deferred`는 DNA 서열 범위의 완료임. `metadata/sample.json`의 `has_rna=true`는 원본 pair 관계이며 `rna_processing_status=deferred`와 구분함. RNA 전체 완료로 집계하지 않음.
- QC 통과 후 이번 실행의 중간 FASTA/개인 VCF는 기존 정책으로 정리함. 이번에 사용하지 않은 기존 WGS/RNA 원본은 보존함. RNA를 나중에 추가할 때 동일 reference/VCF와 저장된 window·좌표를 대조하고 read 근거와 RNA 산출물을 별도로 검증해야 함.
- 정상 RNA pair 처리에서는 보류 옵션을 `0`으로 하고 기존 `MFM_SEQUENCE_RNA_WGS_EVIDENCE=temporary`를 사용함. DNA 전용 완료본을 RNA 전체 완료본으로 재사용하지 않음.

## 샘플별 고정 CPU와 처리 버퍼 설정

`cpu_allocation_mode`는 `fixed_sample` 또는 `dynamic`임. 제공 기본·예제 설정은 fixed_sample/샘플16CPU/동시 최대4/정상 CPU풀68을 사용함. 기존 저장된 실행·사용자 config의 미명시 모드는 호환을 위해 기존 해석을 유지함. 고정 모드는 샘플 시작 시 명시한 CPU 수와 실행 가능한 CPU 집합(affinity)을 확보하고 WGS–RNA 처리 종료까지 유지함. 다른 고정 샘플의 CPU 집합과 분리함. 실제 작업은 기존 `run_sample.sh`와 단계별 도구를 사용하며 별도 제어기는 필요하지 않음.

기존 config의 경로·입력·생물학적 정책을 유지하면서 아래 자원 항목을 적용함. 16×4는 사용 예시이며 최적 조합을 확정한 값은 아님.

```json
{
  "cpu_allocation_mode": "fixed_sample",
  "cpu_utilization": "maximum",
  "cpu_pool_limit": 68,
  "sample_cpus": 16,
  "local_jobs": 4,
  "memory_limit_bytes": 1650000000000,
  "processing_buffers": {}
}
```

CLI로 같은 설정을 선택할 수도 있음.

```bash
# config의 data/output/state 경로와 대상 목록을 먼저 설정함.
bash run_pipeline.sh dry-run --lab medai --batch 001 --config my_config.json \
  --executor local --cpu-allocation-mode fixed_sample \
  --cpu-pool-limit 68 --sample-cpus 16 --local-jobs 4

# 확인한 동일 설정으로 실제 실행함.
bash run_pipeline.sh run --lab medai --batch 001 --config my_config.json \
  --executor local --cpu-allocation-mode fixed_sample \
  --cpu-pool-limit 68 --sample-cpus 16 --local-jobs 4
```

| 샘플당 CPU | 전체68에서 가능한 고정 동시수 예 |
|---:|---:|
| 16 | 최대4 |
| 32 | 최대2 |
| 48/64/68 | 최대1 |

- 소유권: 16×4는64CPU를 유지함. 남은4를 실행 중 샘플에 임의로 더하지 않음. `sample_cpus=auto` 또는 명시 동시수×샘플 CPU가 풀을 넘는 설정은 거부함.
- 하위 배분: 각 도구에 16을 중복 지정하는 방식이 아님. HAP/염색체 worker·정렬·해독·sort·압축·부모 비용의 동시 합계가16 안에 들어가도록 배분함. 도구의 스레드 인자는 실행 직전에 정하며 실행 중 강제로 변경하지 않음.
- 실제 사용: 직렬 도구·파일 읽기/쓰기·최대 worker 수에 의해 실제 사용은 배정보다 작을 수 있음. 다른 프로젝트와 커널 CPU quota를 공유하므로 물리 CPU 독점/항상100% 사용을 보장하지 않음. 다른 프로젝트의 설정은 변경하지 않음.
- 작은 예산: Long의 자동 정렬 최소는 기존4 mapping thread와 보조 비용을 따름. 8CPU로 Long을 실행하려면 기존 실험처럼 `MFM_WGS_MIN_ALIGN_THREADS=3` 등의 최소 스레드 조건을 별도 검토해야 함. 고정 모드만으로 최소 도구 비용을 없애지 않음.
- RAM: 샘플 수로 고정 분할하지 않음. 전역 RAM 운영 상한/10% 여유/입력·버퍼 예측/RSS 보호를 공유함. CPU 부모 소유권은 전역에서 한 번만 세고 하위 CPU는 샘플 내부에서 셈. 부모·모든 단계 RAM은 전역에서 합산함.
- 수명: 정상/실패 종료 시 소유권을 반환함. 부모가 종료해도 살아 있는 분리 단계가 있으면 해당 CPU 소유권을 보존함. 고정 CPU 안에서는 다른 샘플의 유휴 CPU를 추가 대여하지 않음.
- 전환: 실행 중 동적 샘플을 강제 전환하지 않음. 새 비교 실행은 별도 `run_id`, `state_root`, `output_root`로 수행함. Slurm에서는 실제 할당보다 큰 sample CPU를 거부함. 여러 노드를 쓰는 경우에는 노드/할당별 CPU 풀 경로를 분리해야 함.
- 셸 직접 실행: `run_sample.sh`/`run_wgs_sample.sh`/`run_rna_sample.sh`도 `MFM_CPU_ALLOCATION_MODE=fixed_sample`, `MFM_CPU_BUDGET`, `MFM_SAMPLE_CPU_BUDGET`, `MFM_CPU_POOL_DIR`, jj_mfm Python을 명시하면 같은 소유권을 사용함. 일반 사용은 `run_pipeline.sh run`을 권장함.

### 배치·버퍼 옵션

`processing_buffers`를 생략하거나 `{}`이면 기존 기본값을 유지함. 새로운 기본값으로 일괄 축소하지 않음. 같은 환경변수를 서로 다른 값으로 중복 설정하면 오류임.

| 항목 | 연결 인자/변수 | 범위 |
|---|---|---|
| `sort_memory_per_thread` | `WGS_SORT_MEM`, `STAGE04_SORT_MEM` / sort `-m` | 스레드당 근사 버퍼임. 전체 RSS 상한이 아님 |
| `minimap_batch_bases` | minimap2 `-K` | query 배치 크기임. 메모리 감소용500M 이하만 허용함. reference 분할 `-I`는 유지함 |
| `windows_per_chunk` | `WINDOWS_PER_CHUNK` | 처리 배치임. 모델 window 길이·선별 기준은 유지함 |
| `trace_prefetch_records` | `MFM_TRACE_PREFETCH_RECORDS` | 독립 source trace 선행 레코드 수임 |
| `trace_prefetch_packed_mib` | `MFM_TRACE_PREFETCH_PACKED_MIB` | 패킹 선행 버퍼 MiB임 |
| `ase_batch_fragments` | `RNA_ALLELIC_BATCH_SIZE` | ASE fragment 처리 배치임 |
| `star_read_buffer_bytes` | STAR `--limitIObufferSize` | 입출력 각각에 적용할 바이트임. RAM 예측에도 합산함 |
| `star_bam_sort_memory_bytes` | STAR `--limitBAMsortRAM` | STAR BAM 정렬 버퍼임 |

버퍼 축소는 임시파일·merge I/O를 늘릴 수 있음. 같은 입력의 read/변이/window/배열 값과 QC가 동일한지 확인하고 경과시간·RSS·I/O로 판단함. N·blacklist·MAPQ·RNA/ASE/CNV 조건을 바꿔 속도를 얻지 않음.

### CPU·RAM 효율 및 실측 이력

- 대기 공정성: 최소 요구량이 RAM 예산에 맞고30초 이상 기다린 내부/단계 요청은 반환되는 CPU를 모아 입장시킴. 실행 중 자원을 강제 회수하지 않으며 RAM 때문에 들어갈 수 없는 큰 요청은 작은 계산을 막지 않음.
- 신규 동시수: `local_jobs=auto`도 계산해 확정한 J를 신규 입장 상한으로 사용함. 기존 초과 실행은 완료까지 보존함. 내장 큐의 `max_active_samples` 명시값은 계속 사용 가능함.
- 서열 전용 VCF: 신규 Stage03은 최대6CPU의 최소 입장량으로 두 염색체 worker를 사용할 수 있게 함. 샘플/실제 할당이 더 작으면 그 한도를 우선함. 이후 검증 구간은 실제 필요량을 반환함.
- worker RAM: 기본 ProcessPool은 실제 제출 수만큼 생성하는 spawn을 사용함. 살아 있는 유휴 worker의 RAM은 join까지 보존함. 모든 작업이 끝난 유휴 자식만 종료한 뒤 RAM을 반환함. 자식 없는 FASTA 대기는512MiB와 실제 부모 RSS×1.1 중 큰 값을 유지하고 계산 전 피크를 재예약함.
- Stage10: CPU minimap2 종료를 확인한 뒤 mapping CPU만 반환함. 계속 쓰는 sort RAM·다른 HAP 자원은 보존함. 후속 index/QC의 CPU는 기존 단계 경계에서 다시 계산함.
- 관측: 실제 수집 완료 시각을 사용함. 잠금 대기로 오래된 표본은 실제 재수집하며 시각만 갱신하지 않음. `tick_metrics`에 잠금 대기/보유와 판정 시 표본 경과시간을 기록함.
- 잠금: 유효한 기존 RAM 관측은 잠금 없이 재사용함. 동일 worker 장부의 반복 쓰기는 줄이되 CPU/RAM·running/대기 요구량 변경은 즉시 기록함. 살아 있는 내장 감시기가 정상 입장을 맡으면 대기 관리자의 중복 전체 할당은 생략함. 감시기 종료/오래된 상태/PID 재사용/고정 하위 풀은 자체 배정으로 복귀함.
- 이력: 공유 풀의 `resource_metrics.jsonl`에 목표10초 간격으로 CPU 배정/실사용·RAM 예약/RSS·RSS 관측 시각·잠금 시간을 기록함. 실제 간격은 I/O·잠금에 따라 달라질 수 있음. `watchdog.json`은 최신 상태임.
- 단계 감사: 기존 `usage.jsonl`에 고정 CPU 모드/샘플 예산/`affinity_cpu_ids`를 기록함. CPU 사용량은 CPU초/실제 경과시간이며 배정 코어 수와 구분함. RSS 표본은 순간 최대치가 아님.
- 공통 VCF 묶음 추출: 소형 비교에서는 CPU 절감과 임시공간 증가를 함께 관측함. 현재 개별 샘플 재개·정리 계약을 유지하며 다인 공유 임시 패널을 운영 기본값으로 추가하지 않음.


### local-jobs와 cpu-pool-limit의 구분

| 인자 | 의미 |
|---|---|
| `--local-jobs N` | 한 로컬 큐에서 동시에 전처리할 샘플 수 요청임. RAM/profile 검토 후 확정값을 실제 큐와 dry-run에서 동일하게 사용함. process/thread 개수 또는 다운로드 동시수가 아님 |
| `--cpu-pool-limit C` | 같은 CPU 풀을 사용하는 전처리의 정상 CPU 배분 예산 상한임. 샘플마다 C를 배정하는 값이 아님. 실제 affinity/Slurm/cgroup의 더 작은 한도를 우선함 |
| `--sample-cpus S` | 고정 모드에서 샘플 하나가 유지할 CPU 수임. 동적 모드에서는 단계별 실제 grant와 구분함 |
| `--download-jobs D` | 선행 다운로드 동시수이며 local-jobs와 별개임 |

- 예: 고정 모드 `local-jobs4`, `sample-cpus16`, `cpu-pool-limit68`은 최대4개 샘플이 각각16CPU를 소유함. 각 샘플 내부의 여러 worker·도구는 그16 안에서 배분함.
- 실행 중 `max_active_samples` 정책으로 상한을 조정할 수 있음. 고정 모드는 확정 local_jobs보다 늘리거나 null로 해제할 수 없으며, 더 작은 값은 허용함. 동적 모드의 명시 정책 override는 유지함.
- 현재 동적 운영의 정상 예약68에 추가 대여 논리16이 별도로 존재함. 논리 합계84를 실제84CPU 사용으로 해석하지 않음. cpu-pool-limit은 커널 quota를 새로 설정하는 옵션이 아니며 현재 공유 cgroup 실효 CPU 상한은68임.
- 여러 큐가 CPU 예산을 공동 관리하려면 같은 `MFM_CPU_POOL_DIR`를 사용해야 함. local-jobs는 큐별 값이며 서버 전체 샘플 수의 단일 상한이 아님.


### 2026-11-10 기본 실행 정책
- 앞으로 제공 기본 설정과 운영 후속 원본 설정은 최적화된 고정 CPU 모드를 사용함: 샘플16CPU·동시 최대4·정상 CPU풀68·RAM운영1.65TB/예약1.485TB임. 더 작은 큐·RAM·실제 할당 조건은 우선함.
- 별도 자원 인자 없이 run_pipeline.sh run을 실행해도 이 기본 config가 적용됨. 동적 모드는 config 또는 --cpu-allocation-mode dynamic으로 명시 선택 가능함.
- 현재 계산과 effective_config/상태 장부는 보존함. 원본 설정을 다시 읽는 새 batch/run부터 고정 기본을 사용함. 실행 중 작업을 재시작하거나 모드를 강제로 변경하지 않음.


### 실행 중 큐의 후속 샘플만 고정 CPU로 전환

- 공용 풀의 `policy.json`에 아래 운영 정책을 명시하면 이후 시작하는 local worker에만 적용함. Slurm·검증된 GPU 비교 예약 예외와 해당 정책이 없는 별도 풀은 기존 설정을 유지함. 별도 제어기 없이 공식 `run_pipeline.sh run` 경로를 사용함.

```json
{"future_sample_resources": {"cpu_allocation_mode": "fixed_sample", "sample_cpus": 16, "local_jobs": 4}}
```

- 이미 실행 중인 계산·다운로드·부모 제어기는 재시작하지 않음. 이전 제어기는 기존 동시수까지 대기 worker를 둘 수 있으나 고정 CPU 소유권·affinity는 공유 풀 안에서 원자적으로 확보함. 68CPU 풀에서16CPU 고정 샘플은 전체 최대4개이며 더 작은 기존 큐 상한은 유지함.
- `state/worker_resources/<sample>.json`과 샘플 상태의 `cpu_allocation_mode`/`sample_cpus`가 실제 새 worker 정책임. 기존 `effective_config.json`은 이전 제어기의 실행 계약을 보존하므로 이 파일 하나만 보고 현재 worker 모드를 판정하지 않음.
- 고정 CPU 대기는 CPU0·부모 RAM512MiB로 기록함. 30초 이상 기다린 최소 CPU 요구는 기존 공정 배분에 등록하여 다른 작업의 반복 확대에 밀리지 않도록 함. 계산 시작 시 최소4GiB와16CPU를 다시 확인하며 실행 중 예약을 강제로 회수하지 않음.
- 개인 VCF의 `--concat-threads 0`은 최종 압축 스레드 자동값임. 염색체 worker가 모두 종료된 뒤 샘플 CPU 예산의 main1CPU를 제외한 수를 사용함(16CPU이면 압축15). 최종 인덱스도 같은 한도를 사용함. 양수 명시값은 CPU 예산 안에서 보존함. GT·변이 필터·서열·최종 검증은 유지함.
- 정렬/sort 비율은 짧은 부분 실측만으로 기본값을 바꾸지 않음. 도구 스레드·경과시간·CPU초와 read/CIGAR/태그 동등성을 함께 확인함.

- 개인 VCF의 패널 선택과 변환은 별도 병렬 구간임. 16CPU에서는 부모 여유를 둔 선택 최대15개(각1CPU), 변환 최대5개(각최대3CPU)를 사용함. 순서와 bcftools 필터 명령은 유지함.
- 단일 서열 샘플은 계산 전에 단계 RAM을 약13.5GiB 예약하고, 선택 join 후 약5.2GiB, 변환 join 후4GiB로 반환함. 이 값은 단계 예약이며 고정 샘플 부모4GiB는 별도임. 다인 함수 호출은 다른 샘플의 실행을 보호하도록 공통 예약을 유지함.
- 2026-11-10 동일16CPU 부분 실측: 실제 패널 chr1–16 각1Mb의 선택·변환은5.39→3.24초, 개인VCF64MB의 최종압축은3.49→2.11초였음(각2회평균). 변이·GT가 일치함. 두 가속률을 곱하거나 전체 샘플 성능으로 해석하지 않음.
- HG00706 chr22 primary10000reads/1HAP의 정렬·CRAM 비교는 mapping10/sort0 평균84.2초, mapping8/sort2 평균101.6초, mapping6/sort4 평균125.1초였음. CIGAR·태그·복원 레코드와 인덱스 임의 접근을 확인했고 기존 자동 배분 기본값을 유지함.


### 2026-11-11 샘플별 16CPU의 단계 상한 정비

- 고정 샘플 모드에서 WGS·RNA·I/O의 단계별 작은 기본 상한을 제거함. 예를 들어 최종 QC/정리의 기존 최대4CPU 대신 샘플 예산16CPU까지 요청함. 동적 모드의 이전 선택 정책은 별도로 유지함.
- 도구 한 개마다16스레드를 중복 부여하지 않음. 동시에 실행하는 HAP·RNA DAG·입력/정렬/압축 프로세스는 합계16CPU 안에서 배분함. 이 합산 배분에 사용하는 `MFM_CPU_POOL_STAGE_MAX`는 임의 단계 상한과 구분함.
- 고정 모드의 RNA STAR·FASTQ 기본 스레드에 실제 stage grant를 전달함. 명령행 재계산과 HAP/DAG 배분은 기존 내장 실행기를 사용함.
- 최종 QC와 정리 후의 서열 전수검사는 독립 염색체별 프로세스로 실행함.16CPU이면 최대15worker+부모1CPU임. HAP별 저장 배열·N 주석·좌표·내용SHA·원본FASTA 대조 기준은 유지함. 새 helper는 `src/integrate_wgs_rna/validation_parallel.py`임.
- 검증 RAM은 부모4GiB+worker별2GiB를 계산 전에 공유 풀에서 재예약함.15worker의 stage 예상 상한은34GiB이며 고정 샘플 소유권의4GiB는 별도 합산함. worker 주소공간 상한4GiB와 단일 압축 스레드를 적용함. RAM이 부족하면 동시수를 줄이거나 대기함.
- HAP 자산의 SHA 검사는 파일별 스레드로 실행하며 원래 순서로 결과를 병합함. 앞선 검증이 반환한 CPU를 재확보한 다음 읽음. RNA 전장 depth의 합계·nonzero·SHA는 염색체별 병렬 검사하고 창 depth는 기존 순서/청크 검사를 유지함.
- 생물학적 처리·window 필터·최종 합격 기준·삭제 전후 검증 순서는 유지함. 서열 생산 모듈의 코드 지문을 바꾸지 않아 기존 서열을 불필요하게 다시 만들지 않음.
- 현재 고정 worker의 이후 stage 실행과 새 worker부터 적용함. 이미 실행 중인 계산 및 구형 동적 Long 실행은 종료까지 보존함. 공유 커널68CPU·RAM1.65TB/예약1.485TB·GLM 불변 정책은 유지함. I/O 대기·도구 직렬 구간·작업 수 감소 및 다른 프로젝트와의 커널 CPU 경합으로 실제 사용량은16CPU보다 작을 수 있음.

- Stage 12 서열 전용 Zarr의 원본 FASTA 대조는 염색체별로 병렬 실행함. 샘플 16CPU에서 부모 1CPU와 최대 15worker를 배정하며 worker당 RAM 2GiB를 예약함. 염색체 전체가 통과한 뒤 최종 trace와 compact 완료를 기록함. 생산 규격 지문을 보존하여 입력 지문이 같고 원본 대조를 통과한 이전 Zarr를 재사용함.

## Long 서열 및 RNA expression 전용 처리 (2026-11-11)

- Long 서열 전용: `pipeline_env.MFM_LONG_READ_OUTPUT_MODE=sequence_only`를 지정함. Short 옵션과 독립이며 공개 phased VCF·`vcf_consensus`가 필요함. WGS 원본 read 다운로드, 정렬, haplotag, 개인 재정렬, readset/depth 저장을 생략함. 개인 HAP FASTA·chain·window/N/blacklist 주석·서열 Zarr·full/source 검증은 유지함.
- Pair expression 전용: `MFM_RNA_OUTPUT_MODE=coverage_only`, `MFM_SEQUENCE_DEFER_RNA=0`을 지정함. WGS 서열 전용과 함께 사용하면 WGS read 없이 RNA BAM만 추가 수집함. HAP1/2별 GTF lift·STAR index·정렬·duplicate 표시·raw 1bp/retained-window coverage·RNA QC를 수행함. 원래의 MAPQ/duplicate/deletion/mate 처리와 수치형은 보존함.
- ASE/CNV: expression 전용에서 WGS 변이 근거·ASE 집계·CNV/RD·ASE QC는 비적용으로 명시함. 0값의 가짜 통계를 만들지 않음. RNA signal은 HAP 좌표의 총 coverage이며 allele-specific assignment 결과와 구분함.
- 완료 검증: `rna_payload=coverage_only`와 `rna_signal_manifest.json`에 두 HAP QC·annotation·원본 preflight SHA를 기록함. 전장 및 window depth/좌표/내용 SHA 검증을 유지하며 최종검증·삭제전 보호·삭제후검증에서도 ASE 파일을 요구하지 않음. 기존 `full` 모드의 ASE 검증은 유지함.
- 운영: `run_pipeline.sh run --selection <목록.tsv> --config <설정.json> --executor local`로 실행함. CPU는 fixed_sample/16, local_jobs=4, 풀68 안에서 실제 예약함. 다른 payload의 완료본을 덮어쓰지 않도록 별도 output_root/state_root/run_id를 사용함.
- RNA 보류 모드와 coverage_only는 동시 사용하지 않음. 기존 DNA-only 완료본에 RNA를 추가할 때 개인 FASTA·chain이 정리되었다면 같은 VCF/reference/window 계약으로 재생성·대조한 뒤 진행해야 함.

- 서열 전용 내부 정리와 기존 입력 보존: `CLEAN_INTERMEDIATES_ON_SUCCESS=1`로 모델 입력만 남기되, `MFM_MODEL_ONLY_CLEANUP=0`과 `remove_owned_raw_inputs=false`를 사용하면 외부의 기존 raw/VCF는 보존함. 기존 full 모드의 필수 조합은 유지함.
- RNA signal 전용은 WGS Stage02 생략 시에도 FASTQ 생성 전에 동일한 `rna_preflight.py`로 RNA quickcheck·header·SM·contig·index·전수 flagstat·GTF 서명을 검사함. 입력 QC를 생략하지 않으며 기존 preflight가 있으면 같은 파일의 서명 검증을 유지함.

출력 디스크 보호: JSON의 `min_output_free_bytes`는 새 local worker 입장 전에 `output_root`가 위치한 저장소의 여유를 검사함. 미설정/0은 기존 동작을 유지함. 입력이 `/tmp` 등 다른 저장소여도 출력 경로를 별도로 검사하며, 진행 중인 계산을 종료하거나 데이터를 삭제하지 않음. RAM 기반 tmpfs는 영구 저장 디스크로 간주하지 않음.

공개 VCF의 REF가 GRCh38 모호염기를 N으로 표현하는 경우, 검토한 샘플의 `sample_pipeline_env`에 `MFM_REFERENCE_IUPAC_POLICY=to_n`을 명시할 수 있음. 기본은 `preserve`임. 원본 참조를 덮어쓰지 않고 입력 references에 파생 참조를 만들며 A/C/G/T/N·좌표·길이는 유지하고 R/Y 등 IUPAC만 N으로 표현함. source/target 서명과 contig SHA를 기록하고 03/09에 같은 참조를 바인딩함. REF 불일치 검사와 변이 필터는 유지함.

## Pileup feature 전수검증 효율화 (2026-11-18)

- 적용 범위: `pileup_feature`의 독립 source 검증 및 저장 배열 검사임. 생산기의 CIGAR 집계·HAP 서열·window·N/blacklist·read 필터·12채널·삽입 투영 정책은 유지함. 기존 완성본은 같은 입력 지문과 전수검증을 통과하면 재사용함.
- 독립 집계: Python 생산기와 별개인 HTSlib pileup의 native C++ 경로를 기본으로 사용함. base/deletion-span은 pileup으로 집계하고, 같은 anchor의 I와 D가 함께 있는 경우를 보존하도록 사건은 별도 CIGAR 경계로 집계함. BAQ·BQ 필터·mate-overlap 축약·8,000 depth 제한을 추가하지 않음. 상한 도달·입력 손상·지원하지 않는 관측은 실패 처리함.
- 선택: `MFM_PILEUP_VALIDATION_BACKEND=auto`가 기본임. native 빌드/라이브러리 로드 불가일 때만 기존 `aligned_pairs` 검증기로 대체함. `native`는 대체 없이 실패하고, `aligned_pairs`는 기존 검증기를 강제함. 계산·입력 오류는 자동 대체로 숨기지 않음. 현재 conda의 HTSlib와 기존 `native_tool` 빌드·잠금·지문 cache를 사용하며 새 Python 의존성은 없음.
- 읽기 통합: 각 window의 저장 pileup 배열을 한 번 읽어 내용 SHA·유한값·mapping/삽입/padding·모든 source 값을 함께 검사함. 저장 순서의 SHA와 검증 전후 입력 지문 검사는 유지함. 양 HAP의 인접 reference block을 공유하며, 전체 염색체 관측 배열은 만들지 않음.
- 자원: 염색체 worker당 native 집계1CPU, 압축1스레드임. 기존 샘플 CPU 예산과 `cpu_results`의 RAM 보호를 사용함. 16CPU이면 최대15worker+부모1CPU, 단계 RAM 계획은 부모4GiB+worker당2GiB임. 이미 실행 중인 작업·다른 프로젝트 설정은 변경하지 않음.
- 생성과 검증의 메모리 수명: 생성용 reference cache·마지막 window 배열을 source 검증 전에 반환함. worker RAM 예약2GiB와 가상 주소공간 상한4GiB를 분리함. 가상 주소공간에는 라이브러리·allocator·매핑도 포함되므로 RAM 예약을 그대로 RLIMIT_AS에 쓰지 않음. 전체 RAM 예약·실측 RSS 보호는 유지함.

### 2026-11-19 생성·통계·RNA CPU 배분 보완

- Pileup 생산은 `MFM_PILEUP_PRODUCTION_BACKEND=auto`가 기본임. 별도 C++ CIGAR walker를 사용하고 빌드·로드 불가에서만 기존 Python 생산기로 대체함. `native`는 대체 없이 실패하며 `python`은 기존 생산기를 강제함. 관측·필터·품질결측·indel·depth cap 없음·Counts.finish는 유지함. source 검증기는 별도 HTSlib pileup으로 집계하며 생산기의 계수 루프를 호출하지 않음.
- 통계 수집은 완료 compact의 manifest·contig 메타데이터 지문이 같을 때 크기 집계를 재사용함. 생성 중이거나 완료 지문이 바뀌면 새로 측정함. 현재 RNA coverage와 TSV 크기는 매번 실제 파일을 집계하므로 ASE QC 추가분을 반영함. 임의의 chunk 변조는 기존 전수 QC의 책임이며 이 통계 cache를 무결성 검사·삭제 승인에 사용하지 않음. cleanup의 원본 전수 집계는 그대로 유지함.
- 고정 샘플 RNA DAG는 짧은 WGS 변이·CNV 근거 계산 완료 뒤 STAR index를 시작함. 16CPU이면 양 HAP가 각각 최대8CPU를 요청하며 메모리 예약·실제 할당을 우선함. 이미 실행한 STAR의 thread를 중간에 늘리는 방식이 아님. 동적 모드·번호 재개·coverage_only 경로·ASE 선행조건은 유지함.
- 비교 시 readset과 pileup의 Stage12 생성·source 검사 시간을 구분함. 과거 동적 CPU 기록과 새 고정16CPU 기록을 같은 조건의 가속률로 표시하지 않음. WGS01–10·RNA가 공통으로 재사용된 시험은 원본부터 전체 처리시간으로 보고하지 않음.
- 기록: 새 trace manifest의 `source_validation`에 contig별 backend·검증초·window 수·검사 배열 바이트·reference block/cache hit 수를 기록함. 재사용 결과의 검증시간이며 생성시간과 구분함. 샘플링으로 검증 범위를 줄이지 않음.
