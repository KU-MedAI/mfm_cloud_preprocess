# 클라우드 전처리 사용 안내

> 2026-11-24T01:10:35 최신 배정: PoC228명, Azure/AWS 본 처리 각1257명(WGS-only), 서버470명(pair221명 우선·Short WGS-only249명 후속). Long WGS261명은 전원 클라우드다. 이전 배정 인원수의 CSV는 백업 후 현재 경로에서 제거한다. 전송 검증 후 회수 모드는 후속 구현이며 실제 전처리는 시작하지 않았다.


> 2026-11-24T00:49:23 분담 명단 갱신: PoC228명 유지, Azure 본 처리1,160명, AWS 본 처리1,165명, 현재 서버659명. 서버는 Short pair188명을 우선 처리할 계획이며 나머지471명은 기존 WGS-only다. 전체3,212명 중복·누락0. 전송 검증 후 로컬 회수 모드는 아직 구현하지 않았고 실제 전처리·데이터 송신은 실행하지 않았다.


기준일: **2026-11-24**. 현재 `config_v2.sh`, `run_wgs_sample.sh`, `run_rna_sample.sh`, `handoff/` 구현을 기준으로 작성함. 이전 날짜의 실험·운영 기록보다 이 안내와 현재 코드를 우선함.

## 1. 출력 모드 선택

WGS 출력과 RNA 출력은 별개 설정임. JSON의 `pipeline_env`에 값을 지정함. `mode=pileup_feature`라는 단일 최상위 인자는 없으며 `acceleration.mode`는 CPU/GPU 선택임.

| WGS 설정값 | 최종 WGS 입력 | WGS-only 개인 HAP 재정렬(Stage10) | WGS read 입력 |
|---|---|---|---|
| `full` | HAP 서열 + readset sequence/메타데이터 + HAP 기준 depth | 적용함 | 필요함 |
| `sequence_only` | HAP 서열 + window 주석·좌표 매핑 | 비적용임 | 불필요함 |
| `pileup_feature` | HAP 서열 + reference 기준 pileup tensor·depth·mask | **비적용임** | 필요함 |

- 선택 키: `MFM_SHORT_READ_OUTPUT_MODE`, `MFM_LONG_READ_OUTPUT_MODE`임. 둘 다 코드 기본값은 `full`이며 pileup은 명시 선택함.
- 기본 pileup은 **reference 정렬 read → feature 집계 → consensus chain으로 HAP1/2 좌표 투영**임. Stage10 개인 HAP CRAM을 feature 계산에 사용하지 않음. Stage01의 reference 정렬·표준화와 Stage06 haplotagging은 필요함.
- `sequence_only`·`pileup_feature`는 `vcf_consensus`를 사용함. 서열 전용 short는 제공 phased Illumina panel, long은 제공 phased public VCF를 사용함. assembly 직접입력과 혼용하지 않음.
- 모드를 바꿀 때 기존 완료본을 덮어쓰지 않도록 별도 `output_root`·`state_root`·`run_id`를 사용함.

| RNA 설정 | 실행 범위 | WGS 근거 처리 |
|---|---|---|
| `MFM_RNA_OUTPUT_MODE=full` | HAP1/2 RNA 1bp signal + HAP 변이 검증·ASE·CNV | sequence/pileup 출력이어도 WGS 근거·개인 HAP CRAM 처리 필요함 |
| `MFM_RNA_OUTPUT_MODE=coverage_only` | HAP1/2 RNA 1bp signal·window signal·QC, ASE/CNV 비적용 | sequence_only에서는 WGS read 불필요, pileup에서는 feature용 reference WGS read 필요함 |
| `MFM_SEQUENCE_DEFER_RNA=1` | DNA sequence/pileup 먼저 완료하고 RNA 보류 | 비full WGS 출력에서만 적용함. pair 관계와 보류 상태를 metadata에 보존함 |

RNA 기본값은 `full`임. `sequence_only` + RNA `full`은 `MFM_SEQUENCE_RNA_WGS_EVIDENCE=temporary`를 명시해야 함. RNA `coverage_only`와 RNA 보류를 동시에 설정하면 오류임. RNA signal은 HAP 좌표의 총 coverage이며 allele별 expression 분리 결과를 뜻하지 않음.

### 실측 전용 HAP 재정렬 옵션

`MFM_PILEUP_BENCHMARK_REALIGN=0`이 기본임. `1`은 비교 실험을 위해 pileup에도 Stage10 비용을 추가하는 옵션이며, feature 입력을 개인 HAP CRAM으로 바꾸지 않음. **일반 WGS-only pileup 처리에는 0을 사용함.** RNA `full`이 필요로 하는 Stage10은 이 실험 옵션과 별개임. 강제 재정렬 실험 시간을 기본 pileup 시간으로 보고하지 않음.

## 2. 실행 환경과 설정 파일

```bash
cd /path/to/cloud_전처리용
# 새 서버/환경 설치가 필요한 경우에 실행함. 기존 환경 확인만 할 때는 verify-* 사용함.
bash ./00_environment.sh
conda activate jj_mfm
bash ./00_environment.sh verify-core
# RNA를 처리하는 경우 추가 확인함.
bash ./00_environment.sh verify-rna
```

- 새 클라우드에서 `run_pipeline.sh run`은 Conda/Miniforge와 jj_mfm 환경을 자동 준비함.00번을 직접 사용할 경우 Conda 선행 설치가 필요함. 기존 jj_mfm은 검증만 수행하고 자동 변경하지 않음. 부족한 패키지는 `bash run_pipeline.sh env install-core`로 environment.yml/requirements.txt를 반영함.
- 인수 없는 00번은 설치·보강을 수행함. `environment.yml` + `requirements.txt`는 `jj_mfm`, `clair3_environment.yml`은 선택 기능용 별도 `clair3_v2` 환경의 근거임. 설치하지 않고 확인하려면 `verify-core`·`verify-rna` 등을 사용함.
- 도구 버전: `TOOL_VERSIONS.tsv`, Clair3는 `TOOL_VERSIONS_CLAIR3.tsv`임. [환경 도구 안내](환경_도구버전.md) 참조함.
- native pileup은 현재 conda의 HTSlib와 C++ compiler로 빌드·캐시함. `native`를 강제할 경우 빌드 실패를 숨기지 않고 중단함.

```bash
cp handoff_config.example.json my_pileup.json
```

`my_pileup.json`의 기존 항목을 유지하면서 다음 항목을 수정·병합함. 경로는 실행 서버에 맞춤. 상대 경로는 JSON 파일의 위치를 기준으로 해석하므로 다른 폴더에 복사할 때 함께 검토함. 이 복사본은 상위 중앙 저장소를 자동 연결하지 않음. JSON에 지정한 최종 경로를 dry-run에서 확인함.

```json
{
  "data_root": "/data/mfm/raw",
  "output_root": "/data/mfm/output/pileup_v1",
  "state_root": "/data/mfm/state/pileup_v1",
  "run_id": "pileup_v1",
  "cpu_allocation_mode": "fixed_sample",
  "cpu_utilization": "maximum",
  "cpu_pool_limit": 68,
  "sample_cpus": 16,
  "local_jobs": 4,
  "elastic_max_jobs": 4,
  "memory_limit_bytes": 1650000000000,
  "remove_owned_raw_inputs": true,
  "cleanup_shared_fasta": false,
  "pipeline_env": {
    "MFM_SHORT_READ_OUTPUT_MODE": "pileup_feature",
    "MFM_LONG_READ_OUTPUT_MODE": "pileup_feature",
    "MFM_PILEUP_BENCHMARK_REALIGN": "0",
    "MFM_RNA_OUTPUT_MODE": "full",
    "MFM_SEQUENCE_DEFER_RNA": "0",
    "MFM_PILEUP_PRODUCTION_BACKEND": "auto",
    "MFM_PILEUP_VALIDATION_BACKEND": "auto",
    "WINDOW_STORAGE_POLICY": "annotated_acgtn_v1",
    "CLEAN_INTERMEDIATES_ON_SUCCESS": "1",
    "MFM_MODEL_ONLY_CLEANUP": "1",
    "MFM_EXTERNAL_CLEANUP_UNOWNED": "0",
    "CLEAN_SHARED_ASSETS_ON_SUCCESS": "0"
  }
}
```


## 3. 내장 큐로 실행

현재 클라우드 분담은 `--target`으로 선택함. azure_poc228명·azure_main1,257명·aws_main1,257명·server470명임. 명단과 전용 설정을 자동 선택하며 `--config`로 경로/자원을 덮어쓸 수 있음.

```bash
bash run_pipeline.sh dry-run --target azure_poc
bash run_pipeline.sh run --target azure_poc --sample HG00308
```

- `--sample`은 해당 target 안에서 선택함. target과 lab/batch/selection을 함께 지정하지 않음.
- 선택한 명단 전체 완료가 기준임. 기존220명/10TB 목표 중단은 target에 적용하지 않으며 자원/디스크 보호는 유지함.
- 전체target과 단일샘플의 상태키를 분리함. 같은 출력의 동일샘플을 전체 큐와 개별 작업에 중복 배정하지 않음.
- 아래 기존 lab/batch 선택 방식도 호환함.

```bash
# 다운로드·설치·전처리 없이 선택/자원/입력 목록/최종 경로를 확인함.
bash run_pipeline.sh dry-run   --selection manifests/리드별/medai/001/long_read.tsv   --config my_pileup.json --executor local

# 확인한 설정으로 실제 실행함.
bash run_pipeline.sh run   --selection manifests/리드별/medai/001/long_read.tsv   --config my_pileup.json --executor local
```

- `--config`: 경로·자원·출력 모드·정리 정책을 담은 JSON임. 예제를 복사해 사용 가능하며 별도 운영 Python을 작성할 필요 없음. target 실행에서 인자 생략 시 해당 target 예제를 읽고, 기존 lab/batch 실행은 handoff_config.json 또는 루트 기본 예제를 읽음.
- `--executor local`: 현재 서버의 내장 큐를 사용함. `slurm`: config의 Slurm 요청을 사용하며 실제 할당을 우선함.
- Slurm 배포 예제는 샘플 작업마다 `slurm.memory=384G`(384GiB)를 요청함. 두 HAP의 STAR 예약·부모 프로세스·10% 여유를 고려한 값이며 실제 RSS 사용량이 아님. 동시 작업4개면 작업별로 요청하며 local RAM 설정은 바꾸지 않음. 더 작은 사용자 설정은 단계별 예약이 들어가는지 먼저 확인함.
- 기존 `--selection`: 한 기관·한 배치의 선택 TSV임. 이 방식은 1~220명 범위이며 여러 배치 전체를 한 TSV로 합쳐 넣지 않음. `--lab medai --batch 001`로 지정할 수도 있음.
- 기본 실행은 다운로드가 끝난 샘플부터 CPU/RAM·동시수·출력 공간을 확인해 투입함. 다음 샘플 다운로드를 계산과 겹치며 상태는 약1초마다 확인함. 기본 다운로드 재시도와 실패 샘플 격리는 내장 기능임.
- 같은 출력/샘플을 별도 큐에서 중복 실행하지 않음. 기존 lab/batch 방식에서는 `model_target_bytes` 도달과 전체 완료를 구분함. target 방식은 전체 명단 완료가 기준임.
- `run_read_group.sh`는 여러 배치 실행용 보조 진입점임. 개별 처리와 local 큐는 같은 `run_pipeline.sh`를 호출함.

## 4. CPU·RAM·저장공간

| 설정 | 역할 |
|---|---|
| `cpu_pool_limit` / `--cpu-pool-limit` | 같은 풀에서 공유하는 전체 전처리 CPU 배분 상한임 |
| `sample_cpus` / `--sample-cpus` | 고정 모드에서 샘플 하나가 확보하는 논리 CPU 예산임 |
| `local_jobs` / `--local-jobs` | 해당 로컬 큐의 동시 계산 샘플 수임. 다운로드 동시수와 다름 |
| `download_jobs`, `prefetch_samples` | 다운로드 동시수와 미리 준비하는 샘플 수임 |
| `memory_limit_bytes` | 전체 운영 RAM 상한 요청임. 더 작은 Slurm/cgroup/물리 할당을 우선하고 10% 여유를 둠 |
| `min_output_free_bytes` | 각 새 샘플 입장 전에 확인하는 출력 파일시스템 최소 여유임. 미설정/0이면 해당 추가 조건 비활성임 |
| `processing_buffers` | sort/STAR/trace/ASE 처리 버퍼 설정임. 생략/빈 객체는 기존 기본값 유지함 |

- original의 자원 예시는 **고정16CPU × 최대4명·공유풀68CPU·RAM1.65TB·예약1.485TB**임. 클라우드 제공 예제는 샘플16CPU를 유지하고 풀/local_jobs는auto이며 실제 물리/cgroup/Slurm RAM 할당을 반영함. 실제 설정은 [클라우드 실행](docs/클라우드_실행.md)을 우선함. CPU affinity는 같은 고정 풀 안에서 분리하며 다른 프로젝트의 물리 CPU 독점을 보장하지 않음.
- 샘플 내부 HAP·염색체 worker·정렬·해독·압축의 동시 CPU 합계가 샘플 예산 안에 들어가도록 계산함. 스레드 수는 각 명령 시작 전에 결정하며 실행 중 도구의 스레드를 강제로 늘리지 않음.
- `maximum`도 고정 샘플의16CPU를68CPU로 자동 확장하지 않음. **동시수만1로 줄이면 여전히16CPU 예산**임. 한 샘플68CPU 실측은 `--local-jobs 1 --sample-cpus 68 --cpu-pool-limit 68`을 함께 지정함. 기존 실행 도중 조건을 섞지 않음.
- 고정 모드는 샘플 CPU 양의 정수를 사용함. `sample_cpus=auto`는 동적 모드의 자동 계산과 구분함. RAM은 샘플 수로 균등 분할하지 않고 단계별 입력·버퍼·worker 피크에 따라 예약·반환함.
- 최근 HiFi16CPU 실측의 자동 배분 예는 HAP1/2 순차·minimap2 -t10·sort -@0·해독 보조1과 기타 파이프라인 비용임. 각 도구에16씩 주지 않으며 HAP 동시수·스레드는 CPU/RAM·최소 정렬 스레드 조건으로 계산함.
- `samtools sort -m`은 스레드당 근사 버퍼임. `-@ 0`도 주 스레드로 실행함. 버퍼를 줄이면 임시파일·merge I/O가 늘 수 있으므로 전체 RSS/시간/I/O를 함께 확인함.
- pileup worker 기본값은 RAM 예약4GiB와 가상 주소공간 보호8GiB임. `MFM_PILEUP_WORKER_RAM_BYTES`, `MFM_PILEUP_WORKER_AS_BYTES`로 명시 가능하며 주소공간은 RAM 예약 이상이어야 함. 가상 주소공간·예약·실제 RSS는 다른 지표임.
- 다수 큐가 같은 노드의 예산을 공유하려면 같은 `MFM_CPU_POOL_DIR`을 사용함. 서로 다른 노드/Slurm 할당의 CPU 풀은 분리함.

## 5. 모드별 실행 단계

| 단계 | 기능 | 모드에 따른 차이 |
|---|---|---|
| WGS01–02 | canonical BAM 준비·입력 검증 | WGS read 불필요한 sequence_only에서는 비적용임 |
| WGS03 | sequence/tag VCF 준비 | 제공된 개인 VCF 또는 Illumina panel을 사용함 |
| WGS04–05 | Clair3 self-call/QC | 명시적 opt-in일 때만 적용함 |
| WGS06 | reference read haplotag | full/pileup 및 RNA full의 임시 WGS 근거에 적용함 |
| WGS07–08 | de novo/anchor | 해당 모드 opt-in일 때만 적용함 |
| WGS09 | HAP1/2 FASTA·FAI·chain 생성 | sequence/pileup도 유지함 |
| WGS10 | 개인 HAP WGS 재정렬·CRAM 출력 | full, RNA full용 근거, 강제 실험 옵션에서만 적용함 |
| WGS11 | window·N/blacklist 주석 | HAP CRAM 검사는 Stage10이 필요한 경로에서만 수행함 |
| WGS12 | 최종 Zarr·전수 source 검증 | full=readset, sequence_only=DNA, pileup_feature=DNA+reference pileup임 |
| WGS23–24 | 최종 WGS QC·선택적 정리 | WGS-only 또는 RNA 보류에 적용함. pair 실행 중에는 RNA 완료까지 미룸 |
| RNA13–19 | FASTQ·lift·STAR·중복표시·1bp signal·QC | RNA full/coverage_only 모두 적용함 |
| RNA20–23 | WGS 변이 근거·ASE·CNV·ASE QC | RNA full에 적용함. coverage_only는 비적용임 |
| RNA24–26 | metadata·최종 QC·선택적 정리 | RNA payload에 맞는 필수 산출물을 검증함 |

대표 도구와 각 단계의 입력/출력은 [PIPELINE](PIPELINE.md), [단계별 상세](docs/단계별_입출력_상세.md)에 있음. RNA full은 의존성 DAG와 HAP 병렬화를 사용하며, 단계 재개·중단값을 명시한 경우 실제 로그의 실행 순서를 확인함.

## 6. 서열·window와 pileup 해석

- GRCh38 + phased SNV/indel VCF로 HAP1/2 생성함. 모델 염색체는 chr1–22임.
- Window는 1,048,576bp 비중첩이며 seed42의 결정적 시작점을 사용함. HAP/염색체 양끝의 연속 N만 trim하고 내부 N은 보존함. 실제 구간에 ACGT가 없는 창은 제외하며 padding은 구분함.
- Blacklist로 창을 삭제하지 않음. blacklist 위치·겹침 비율과 N 수·비율·연속 길이 등을 window metadata에 기록함. N 1%/5% 적격성 같은 중복 boolean 대신 `n_fraction`으로 학습 시 기준을 선택함.
- Pileup은 HP별/allele별 지지 tensor가 아님. 같은 reference 위치에 대응하는 HAP1/2는 같은 관측을 사용함. HAP 서열은 각 allele을 보존함.
- HAP insertion은 서열에 그대로 있고 `mapping_kind=2`임. 왼쪽 reference anchor의 generic insertion 사건 지지를 전달함. 삽입 allele 일치율이나 삽입 염기별 실제 depth/BQ/base ratio로 해석하지 않음.
- DNA 복원 정답에서는 N/padding P를 제외하고 유효 정답 수로 정규화함. RNA/WGS 실제 관측값은 N·blacklist 위치라는 이유만으로0으로 덮어쓰지 않음. 이 계약은 window metadata에 기록하며 학습 코드가 적용함.
- 저장: HAP DNA 4bit ACGTNP, pileup 12채널 float32, count/depth uint32, zstd level5임. MAPQ·primary/supplementary·duplicate 정책은 플랫폼의 기존 필터를 따름. 상세 배열/채널은 [pileup 상세 설명](docs/전처리_상세.md#reference-pileup-선택-모드-2026-11-23)에 있음.

## 7. 검증·재개·정리

- full compact 검증 + 독립 source trace + 최종 패키지 QC를 유지함. native 생산은 C++ CIGAR walker, native 독립 검증은 별도 HTSlib pileup 집계임. 검증 범위를 줄여 속도를 얻지 않음.
- Backend `auto`는 빌드/로드 불가에서만 Python/기존 aligned-pairs로 대체함. 입력 손상·집계 오류·불일치는 실패하며 fallback으로 숨기지 않음.
- 입력/정책 지문이 같은 완료 contig만 전수 대조 후 재사용함. 다른 payload나 지문으로 기존 Zarr를 덮어쓰지 않음.

```bash
# 실패 원인을 해결하고 선행 자료가 남아 있는지 확인한 뒤 요청함.
# 기존 선택 TSV와 config를 그대로 사용함. 큐가 종료돼 있으면 같은 run 명령을 다시 실행함.
bash run_pipeline.sh resume   --selection manifests/리드별/medai/001/long_read.tsv   --config my_pileup.json --sample NA12878 --from-stage 12
```

내장 WGS 재개 요청은 검토된 단계01/02/06/10/12/23을 지원함. 직접 branch 실행의 `WGS_START_AT`·`WGS_STOP_AFTER`, `RNA_START_AT`·`RNA_STOP_AFTER`는 해당 branch 번호를 사용함. RNA 번호를 WGS `--from-stage`에 넣지 않음.

| 보존 목적 | 설정·동작 |
|---|---|
| 최종 모델 입력 중심 보존 | `CLEAN_INTERMEDIATES_ON_SUCCESS=1`, `MFM_MODEL_ONLY_CLEANUP=1` 사용함. 최종 QC와 삭제 전 guard 통과 뒤 내부 중간 정리·정리 후 검증을 수행함 |
| 기존 외부 원본은 남기고 새 출력 내부 중간만 정리 | 비full 출력에서 `CLEAN_INTERMEDIATES_ON_SUCCESS=1`, `MFM_MODEL_ONLY_CLEANUP=0`, `remove_owned_raw_inputs=false` 사용함. 외부 cleanup을 차단함 |
| 재실측용 중간자료 보존 | `CLEAN_INTERMEDIATES_ON_SUCCESS=0`, `MFM_MODEL_ONLY_CLEANUP=0`, `remove_owned_raw_inputs=false` 사용함. 조각/run BAM까지 필요하면 `MFM_HAPLOTAG_KEEP_SHARDS=1`, `STAGE04_KEEP_RUN_BAMS=1`도 명시함 |

**`remove_owned_raw_inputs=false`만으로 외부 정리 전체를 차단하지 않음.** `MFM_MODEL_ONLY_CLEANUP=1`은 환경이 갖춰지면 별도 `completion_cleanup`을 호출하므로, 기존 외부 raw를 반드시 보존할 실행은 위 조합을 함께 적용함.

모델 전용 정리는 `compact/`, 최종 RNA Zarr/TSV, `metadata/`, `qc/`, `logs/`를 보존함. FASTA·비압축 chain·개인 BAM/CRAM·STAR index·FASTQ·작업폴더는 정리 가능함. window 원장과 압축 chain은 `metadata/wgs_windows*`, `metadata/reference_mapping/`에 남김. FASTA 파일을 항상 영구 보존한다고 가정하지 않음. 공유 reference/panel·다른 샘플·소유권 미등록 원본은 자동 일괄 삭제하지 않음.

최종 QC가 없거나 FAIL/stale이면 정리하지 않음. RNA full은 ASE/CNV 근거도 확인하고 coverage_only는 해당 결과를 요구하지 않음. 중간 정리 뒤 RNA를 추가/재실행하려면 같은 reference/VCF로 개인 자산·필요한 read 근거를 재구성하고 window/서열 일치를 검증해야 함.

## 8. 최근 효율화와 구현 위치

| 구현 | 코드 위치·유지하는 계약 |
|---|---|
| 고정 샘플 CPU·단계별 재배분 | `src/integrate_wgs_rna/sample_cpu.py`, `scripts/resource_runner.py`, `lib/stage_runner.sh`임. 소유CPU·도구 합산·메모리 보호를 유지함 |
| HAP 재정렬 CPU/RAM 계산 | `src/integrate_wgs_rna/realign/align_resources.py`, `align_memory.py`, `phases.py`임. 실행 중 정렬기 강제 변경 없이 경계에서 자원을 재확보함 |
| Pileup 생성·독립 검증 | `pileup_features/compute.py`, `native_oracle.py`, `projection.py`, `export.py`임. 필터·12채널·삽입 투영·전수 source 대조 유지함 |
| Pileup worker 메모리 보호 | `export.py::_worker_memory_limits`임. RSS와 주소공간을 구분해 기본예약4GiB/주소8GiB 사용함 |
| 완료 compact 크기 집계 | `handoff/metrics.py`임. 완료 manifest/contig 지문이 같을 때 집계를 재사용하며 전수 QC나 삭제 승인으로 대체하지 않음 |
| RNA 병렬화·입력/검증 재사용 | `lib/rna_dag.sh`, `rna_cache.py`, `rna_profile.py`임. 실제 grant·입력/설정 지문·ASE/CNV 의존 관계를 유지함 |

native backend/버퍼·worker 동시수를 바꿔도 서열·window·read 필터·수치·QC 범위가 같아야 함. 작은 구간의 가속률과 전체 샘플 시간을 구분하며 생물학적 기준을 완화해 속도를 얻지 않음.

## 9. 출력과 상태 확인

- 출력 경로: [실행·출력 구조](docs/실행과_출력구조.md) 참조함. `run_root`는 지정 root 바로 아래, 기본 배치 layout은 그 아래 `<lab>_<batch>`를 붙임.
- 실제 확정 config: `<state_root>/<lab>_<batch>/effective_config.json`임. 원래 config만 보고 실행 중 worker의 조건을 단정하지 않음.
- 큐 상태: 같은 상태 폴더의 `queue/status.json`, 샘플 상태는 `samples/<sample>.json`·`samples/<sample>.run.log`임. 단계별 상세 provenance/QC는 샘플의 `metadata/`·`qc/`·`logs/`에 있음.
- 자원 보고는 배정 CPU, 실제 CPU초/경과초, 관측 RSS/PSS, 프로세스 I/O를 구분함. `/proc/PID/io`는 Lustre 물리 디스크 전체 속도가 아님. 10초 관측최대는 순간최대 보장이 아님.
- 실측 시 다운로드·수동정지·입력 재사용·재시도·정리 여부를 함께 기록함. Stage12 부분 시간과 샘플 전체 시간을 혼용하지 않음. readset과 pileup의 CPU 배분·동시수·필수 단계가 다르면 저장 형식만의 가속률로 해석하지 않음.
