# MFM 클라우드 전처리

**`cloud_전처리용`의 전처리 코드·환경 정의·입력 목록·Docker 구성을 함께 포함한 저장소임.** 이 폴더 전체로 빌드하며, 빌드 중 기존 비공개 소스 저장소를 추가로 내려받지 않음. 전처리 알고리즘·단계 번호는 최신 `original` 기반을 유지함.

- Docker 사용: [Docker 실행](docs/Docker_실행.md)에 따라 이미지 빌드 → 저장 경로·JSON 연결 → `run`을 실행함. VM에 Conda와 분석 도구를 따로 설치할 필요가 없음.
- Docker 없이 사용: 아래 `run_pipeline.sh run`이 Miniforge·`jj_mfm`과 필수 도구를 자동 준비함. 기존 환경이 있으면 검증하며 자동 업데이트하지 않음.
- 입력 준비: 공개 reference·VCF·유전자 주석·샘플 원본은 실행 중 필요한 항목을 다운로드하거나 사전에 준비한 파일을 재사용함. 네트워크·CPU/RAM·충분한 작업 디스크가 필요함.
- 검증 범위: 코드·별도 Conda 환경·작은 실제 입력의 회귀검사를 수행함. **실제 Docker 이미지 빌드와 새 VM의 대용량 샘플 전과정 완주는 아직 검증하지 않음.**
- 현재 배정: PoC228명, Azure/AWS 본 처리 각1,257명, 서버470명으로 총3,212명임. VM 생성·결과 자동 반출은 이 코드의 전처리 실행과 별도임.

## 1. 먼저 확인할 문서

| 문서 | 내용 |
|---|---|
| [클라우드 실행](docs/클라우드_실행.md) | 환경 설치, VM/노드 경로, 단일 샘플/배치, 결과 보존·반출·재개 |
| [Docker 실행](docs/Docker_실행.md) | CPU 이미지 빌드·외부 데이터·CPU/RAM·중단/재개·검증 한계 |
| [pileup·readset 비교](docs/pileup_필터_비교.md) | 공통 read 필터·Short/Long 차이·정렬 좌표·삽입 처리·RNA 조건 |
| [GUIDE](GUIDE.md) | 출력 모드, 자원, 검증/정리 설정과 구현 위치 |
| [PIPELINE](PIPELINE.md) | original 번호 그대로인 WGS01–12/23–24·RNA13–26 단계 |
| [출력 구조](docs/실행과_출력구조.md) | Zarr·window 주석·압축 좌표 매핑·QC·상태 장부 |
| [SOURCE_MANIFEST.json](SOURCE_MANIFEST.json) | 복사한 original 파일의 SHA256와 제외 범위 |
| [검증 결과](docs/검증결과.md) | 회귀·경로/할당·모드·재포장 검사와 실제 클라우드 검증 한계 |

- 코드·목록·작은 blacklist 주석을 포함함. 원본 BAM/CRAM·완료 Zarr·실행 장부·캐시·설치된 GPU 바이너리·수십GB 공통 자산은 포함하지 않음.
- 공통 FASTA/VCF/GENCODE는 동봉 [공통 자산 목록](manifests/공통자산.json)의 공개 URL·파생 규칙·SHA256로 필요한 항목을 자동 준비함. 이미 검증된 `bundled_assets_root`를 재사용할 수도 있음. 목록57개는 메타데이터이며 대용량 실물 파일은 별도로 수집함.
- 기존 연구실 데이터·`original` 실행과 연결되는 외부 symlink나 상위 저장소 자동 탐색을 사용하지 않음. 입력·출력·상태 경로를 JSON에 명시함.
- Azure VM/Batch task 등에서 Linux 명령으로 실행할 수 있는 코드 구성임. Azure 리소스 생성·업로드·실제 노드 전량 실행은 별도 작업임.

## 2. 빠른 시작

```bash
cd /path/to/mfm_cloud_preprocess
# 코드 파일 검사에는 Python 3.10 이상만 필요함.
python3 tools/check_cloud_package.py

# full/readset을 사용할 때는 handoff_config.example.json을 복사함.
# 이번 reference pileup을 사용할 경우 아래 예제를 복사함.
cp cloud_config.pileup.example.json my_config.json
```

`my_config.json`의 data_root/output_root/state_root를 충분한 POSIX 디스크 경로로 수정함. 상대 경로는 설정 파일 위치 기준임. 샘플 기본16CPU와 노드에 맞는 RAM 예산을 확인함.

검증 대상은 Linux x86_64임. Bash·Python 3.10 이상·curl·sha256sum 등 기본 OS 도구, 외부 HTTPS/FTP 접속, 쓰기/실행 가능한 디스크가 필요함. **아래 `run`은 `jj_mfm`이 없으면 Conda/Miniforge와 환경을 자동 준비함.** 기본 Miniforge 설치 경로는 코드 폴더의 `.conda_miniforge`이며 `MFM_CONDA_ROOT`로 변경 가능함. 연구실의 Conda 설치 경로는 사용하지 않음. 최초 환경 설치와 원본/공통 입력 다운로드에는 시간이 필요함.

환경 자동화 구현은 [handoff/environment.sh](handoff/environment.sh)이며 [run_pipeline.sh](run_pipeline.sh)의 run/prepare/worker/finish에서 호출함. **기존 jj_mfm은 검증하며 자동 업데이트하지 않음.** 부족한 환경은 `bash run_pipeline.sh env install-core`로 명시 보강함.

`doctor`·`dry-run`·`download --dry-run`은 전처리 패키지 설치 전에도 확인 가능함. Conda가 이미 설치되어 있고 환경을 미리 설치/보강하려면 `bash 00_environment.sh install-core`를 사용함. 인수 없는00번은 core와 선택 Clair3 환경을 함께 설치하며,00번 자체는 Conda가 먼저 필요함. 기본 제공3212명 목록의 public/panel VCF 처리는 Clair3 self-call을 사용하지 않음.

```bash
# Azure PoC 전체 명단228명과 전용 설정을 자동 선택함. 실제 계산은 하지 않음.
bash run_pipeline.sh dry-run --target azure_poc
# VM 한 대에서 한 명을 처리하는 예시임. 사용자 경로/자원 설정을 적용함.
bash run_pipeline.sh dry-run --target azure_poc --sample HG00308 --config my_config.json
bash run_pipeline.sh run --target azure_poc --sample HG00308 --config my_config.json
```

- `--target`: azure_poc228명·azure_main1,257명·aws_main1,257명·server470명 중 선택함. `--config`가 없으면 해당 target의 설정을 자동 사용함. `--sample ID`로 명단 안의 샘플만 선택 가능함.
- target 실행은 선택한 명단 전체 완료를 기준으로 하며 기존220명/10TB 중단 기준을 적용하지 않음. CPU/RAM/디스크 보호는 유지함. 전체 명단 선택이 전체 동시 실행을 뜻하지 않음.
- target과 `--lab/--batch/--selection`은 함께 사용하지 않음. 기존 lab/batch/selection 방식도 그대로 지원함.
- `--lab`/`--batch`는 동결된 샘플 배정 표의 식별자이며 실행 VM의 소유 기관을 뜻하지 않음.
- `--selection`은 한 기관·한 배치의 TSV 전체를 선택함. **run에서 `--selection`과 `--sample`을 함께 지정해 한 명으로 줄이려 하지 않음.** 한 명은 위의 lab/batch/sample 방식 또는 한 행의 selection을 사용함.
- 다운로드가 끝난 샘플부터 CPU/RAM·디스크 조건 안에 즉시 투입하며 다음 다운로드도 이어가는 내장 큐임. 별도 운영 제어 Python은 필요하지 않음.
- `run_handoff.sh`는 호환 alias임. 직접 branch/과거 목록/Slurm/smoke 명령은 `scripts/run_compat_pipeline.sh`에 보존함. 기본 사용은 `run_pipeline.sh dry-run/run/resume`임.

## 3. 출력 모드와 자원 기본값

아래 표는 루트 예제 기준임. `--target`의 자동 설정은 분담/설정 폴더의 pileup/RNA full 정책이며, 클라우드3예제는 샘플16CPU·동시1/중간보존, server 예제는 공유68CPU·최대4/중간보존임. 설정파일을 명시하면 해당 값을 우선함.

| 항목 | 클라우드 예제 설정 |
|---|---|
| `handoff_config.example.json` | 원래 `full` WGS readset 및 RNA full 정책 유지함 |
| `cloud_config.pileup.example.json` | WGS pileup_feature, RNA full, 강제 HAP 재정렬 옵션0임 |
| CPU | fixed_sample/샘플16, 풀·로컬 동시수 auto임. 노드의 affinity/cgroup/Slurm 할당을 우선함 |
| RAM | 예제 요청1.65TB이되 실제 물리/cgroup/Slurm 한도로 낮춰 반영하고10%여유를 둠. 사용자가 더 작은 운영 예산을 명시 가능함 |
| 정리 | 최종QC/독립trace 후 내부·소유 원본 정리와 정리후 검증함. 기존 외부 raw 보존은 GUIDE의 조합을 사용함 |
| GPU | 기본CPU임. 기존 선택 backend 코드를 포함하되 컨테이너/드라이버/호환 도구는 노드에 별도 설치·검증함 |

- 샘플 예산보다 VM의 CPU가 작으면 고정 모드는 실패함. 샘플CPU/동시수를 config 또는 CLI로 명시 조정함. 단독68CPU는 `--local-jobs 1 --sample-cpus 68 --cpu-pool-limit 68`을 함께 사용함.
- 풀 경로는 같은 노드/할당의 큐만 공유함. 다른 노드의 PID·CPU를 같은 풀에서 예약하지 않음. 각 Batch task는 별도 state_root/run_id를 사용함.
- `MFM_SCRATCH_ROOT`는 지원 단계 임시파일용임. 모든 BAM·STAR index·최종 출력을 NVMe로 자동 이동하는 설정은 아님. 전체 임시작업을 로컬에서 처리하려면 JSON 경로도 로컬 디스크에 지정함.

## Reference pileup 선택 모드 (2026-11-23)

**WGS-only/coverage_only의 기본 pileup에는 개인 HAP WGS 재정렬이 필요하지 않음.** Stage06 reference haplotag BAM에서 집계한 뒤 chain으로 개인 HAP window에 투영함. `MFM_PILEUP_BENCHMARK_REALIGN=0`을 사용함. RNA full의 ASE/CNV가 요구하는 WGS 근거 처리는 별개임.

저장 위치: `compact/hap_sequence_reference_pileup_v1/<chr>/hap_sequence.zarr` 및 `pileup.zarr/HAP1|HAP2`임. 기존 서열·window·N/blacklist·padding·순서·좌표 주석을 유지하고 WGS readset 대신 feature를 저장함.

| 채널 | 의미 |
|---|---|
| `log1p_depth` | primary 염기 관측 depth의 log(1+depth)임 |
| `mean_mapq`, `mean_baseq` | 염기 관측의 평균 MAPQ/BQ임. BQ 결측은 baseq_depth로 구분함 |
| `strand_forward_ratio` | 관측 중 forward 비율임 |
| `read_A/C/G/T/N_ratio` | 실제 base depth에 대한 각 염기 비율5개임 |
| `insertion_ratio` | 왼쪽 anchor 뒤 insertion 사건수/(anchor base depth+deletion span depth)임 |
| `deletion_after_event_ratio` | anchor 뒤 deletion 사건수/동일 분모임 |
| `deletion_after_max_len` | anchor 뒤 read CIGAR deletion의 최대 길이임 |

12채널 float32 외에 uint32 depth/baseq_depth/supplementary_depth/max_deletion_length, uint8 mapping_kind, offsets/lengths를 저장함. HAP 삽입 서열은 DNA에 남기며 `mapping_kind=2`에 왼쪽 anchor의 generic insertion 지지만 전달함. 해당 삽입 allele 일치율/염기별 실제 BQ/depth는 생성하지 않음. 같은 reference 위치의 양 HAP는 같은 관측을 쓰므로 allele별 tensor로 해석하지 않음.

기본 생산기 C++ CIGAR walker와 독립 HTSlib source 검증을 유지함. worker RAM예약4GiB·주소공간 보호8GiB이며 RAM 예산에 맞춰 병렬도를 조정함. 생성·검증 범위를 줄이지 않음.

## 4. 코드 점검과 재배포

```bash
# 데이터 전송/설치 없이 코드와 배포 파일의 해시를 확인함.
python tools/check_cloud_package.py
# 새 경로에 코드·목록·문서만 내보냄. 실행 후 쌓인 data/output/run_state는 포함하지 않음.
bash scripts/make_handoff_bundle.sh /new/path/cloud_bundle
```

`CLOUD_MANIFEST.json`은 배포 시점의 코드/문서/목록 inventory임. 설정은 예제를 복사한 새 파일에서 변경함. 원본 코드나 예제를 바꾸면 inventory를 갱신한 새 배포본을 만들고 변경 근거를 남겨야 함. 이 해시는 전처리 결과 QC나 Azure 실행 성공을 대신하지 않음.

과거 원본의 상세 연구/성능 이력은 [original 변경 이력](docs/original_변경이력.md)에 보존함. 그 문서의 서버 경로·옛 RAM 정책·실험 옵션을 클라우드 기본값으로 사용하지 않음.
