# Docker로 전처리 실행

## 무엇이 달라지는가

- 설치: 코드와 jj_mfm·STAR/samtools·C++ 개발 도구를 하나의 CPU 이미지로 묶음. 서버마다 Conda 환경을 다시 구성하지 않아도 됨.
- 재현: 동일한 이미지 digest를 사용하면 같은 실행 파일·라이브러리를 배포 가능함. Dockerfile 재빌드 시에는 고정하지 않은 전이 의존성과 OS 패키지가 갱신될 수 있으므로 새 이미지를 다시 검증함.
- 배포: 같은 이미지를 Azure/AWS/연구실 Linux 노드에서 실행 가능함. 데이터 위치·CPU/RAM·동시수는 JSON과 Docker 실행 인자로 지정함.
- 속도: Docker 자체가 정렬·pileup 계산을 가속하는 것은 아님. 실제 CPU/GPU·RAM·스토리지와 내부 병렬화가 처리시간을 결정함.
- 결과: 현재 original 기반 생물학 알고리즘·window·pileup·RNA·검증 정책을 유지함. Docker 구성은 실행 환경과 데이터 연결을 담당함.

## 1. 이미지 만들기

Docker Engine이 실행되는 Linux x86_64 호스트에서 이 코드 폴더를 사용함. 최초 빌드는 Ubuntu 저장소·GitHub·Conda 채널 접속이 필요함.

Dockerfile·전용 JSON·entrypoint·소형검사는 `docker/`에서 관리함. 공통 환경 정의와 `.dockerignore`는 루트에 두며, 빌드할 때 마지막 `.`으로 저장소 전체를 context로 지정함.

```bash
cd /path/to/mfm_cloud_preprocess
python3 tools/check_cloud_package.py
# 저장소 루트를 build context로 사용함. Docker 전용 파일은 docker/에 있음.
docker build --platform linux/amd64 -f docker/Dockerfile -t mfm-preprocess:261124-cpu .
docker run --rm mfm-preprocess:261124-cpu
```

- Ubuntu24.04 베이스는 digest, Miniforge26.7.2-0 installer는 SHA256를 고정함. 환경 정의 레이어와 코드 레이어를 나누어 코드 수정만으로 환경 전체를 재설치하는 일을 줄임.
- Python 패키지는 environment.yml/requirements.txt를 사용함. 빌드 중 pip check·코드 해시·00번 core/RNA 검사·native13종 빌드·소형 WGS/RNA와 자원 설정13검사를 수행하며 실패 시 build도 실패함.
- bare run은 도움말만 표시함. 실제 작업은 아래처럼 명시한 run 인자로 시작함.
- 기본 UID:GID는10001:10001임. 호스트 결과 소유자를 맞추려면 아래 --user 옵션을 사용함.
- 기존 클라우드 코드 폴더에 별도 운영파일/자격정보를 추가했다면 먼저 `python3 tools/check_cloud_package.py --export-to /path/to/new-code-only-context`로 코드만 내보내고 그 경로에서 build함. 자격정보나 원본 자료를 이미지에 추가하지 않음.

## 2. 데이터와 설정 연결

| 컨테이너 경로 | 용도 | 보존 |
|---|---|---|
| /opt/mfm | 코드·manifest | 이미지에 포함함 |
| /opt/conda/envs/jj_mfm | 실행 환경 | 이미지에 포함함 |
| /data/input | 개인 입력·수신 파일 | 호스트에 보존함 |
| /data/assets | 공통 reference/panel/GENCODE | 호스트 공통캐시로 재사용함 |
| /data/output | 최종 모델 입력·QC·로그·현재 작업자료 | 호스트에 보존함 |
| /data/state | 실행 상태·재개 원장 | output과 함께 보존함 |
| /scratch | 임시폴더·현재 코드가 지원하는 단계별 scratch | 호스트 로컬 디스크를 사용함 |
| /cache | native 바이너리 캐시 | 호스트에서 재사용 가능함 |
| /config/job.json | 실행할 JSON 설정 | 읽기 전용으로 연결함 |

아래 명령은 Bash에서 실행함. 데이터 폴더는 빌드 코드 폴더 밖에 생성함.

```bash
MFM_DOCKER_ROOT="$PWD/../mfm_docker_run"
mkdir -p "$MFM_DOCKER_ROOT"/data/{input,assets,output,state}
mkdir -p "$MFM_DOCKER_ROOT"/scratch "$MFM_DOCKER_ROOT"/cache
cp docker/config.pileup.example.json "$MFM_DOCKER_ROOT/config.json"

# 여러 docker run에서 재사용하는 실행 옵션임.
MFM_DOCKER_ARGS=(
  --rm --user "$(id -u):$(id -g)"
  --cpus 16 --memory 384g --stop-timeout 120
  --mount "type=bind,src=$MFM_DOCKER_ROOT/data,dst=/data"
  --mount "type=bind,src=$MFM_DOCKER_ROOT/scratch,dst=/scratch"
  --mount "type=bind,src=$MFM_DOCKER_ROOT/cache,dst=/cache"
  --mount "type=bind,src=$MFM_DOCKER_ROOT/config.json,dst=/config/job.json,readonly"
)

docker run "${MFM_DOCKER_ARGS[@]}" mfm-preprocess:261124-cpu \
  doctor --config /config/job.json

docker run "${MFM_DOCKER_ARGS[@]}" mfm-preprocess:261124-cpu \
  dry-run --target azure_poc --sample HG00308 --config /config/job.json

docker run "${MFM_DOCKER_ARGS[@]}" --name mfm-hg00308 mfm-preprocess:261124-cpu \
  run --target azure_poc --sample HG00308 --config /config/job.json
```

- Docker에서는 `--config /config/job.json`을 명시함. 내장 target JSON의 기본 runtime 경로는 일반 설치용이며 컨테이너의 /data 마운트 예제와 다름.
- 294명/전체3212명 등의 명단 배정은 기존 --target/--selection/--lab/--batch를 그대로 사용함. 위 명령의 --sample을 빼면 해당 target 전체를 내장 큐가 순차/병렬 처리함.
- 제공 Docker JSON은 reference pileup, RNA full, full QC/source trace, 최종 QC 뒤 모델전용 정리 정책임. 소유한 다운로드 원본·내부 중간을 정리하며 기존 외부 원본의 무조건 삭제는 허용하지 않음. 보존이 필요한 실험은 기존 GUIDE의 cleanup 설정을 바꿈.
- cache 경로는 바이너리 실행·공유 라이브러리 매핑과 flock을 지원해야 함. noexec 임시폴더나 객체 저장소 경로를 native cache로 사용하지 않음.
- 공통 입력이 없으면 기존 downloader가 준비하므로 input/assets에도 쓰기 권한이 필요함. 이미 검증된 공통파일을 같은 마운트 경로에 준비하면 재사용함.
- /scratch가 모든 중간자료의 저장 위치를 강제로 바꾸는 것은 아님. 현재 scratch 옵션이 적용되지 않은 RNA STAR 작업자료 등은 output 아래 _work를 사용할 수 있음. output 호스트 디스크도 중간자료의 피크를 수용해야 함.

## 3. CPU·RAM 설정

- Docker --cpus/--memory: 컨테이너가 사용할 수 있는 커널 한도임.
- JSON cpu_pool_limit: 내부 전처리 CPU 풀 상한임.
- JSON sample_cpus/local_jobs: 고정 샘플 CPU와 동시 샘플 수임.
- JSON memory_limit_bytes: 내부 RAM 운영 상한임. 실제 cgroup/물리 한도가 더 작으면 우선하며10%여유를 유지함.
- 예제는 CPU16·동시1·RAM384GiB임.384GiB는 두 HAP STAR 예약·부모·여유를 고려한 검증용 설정이며 실제 RSS를 뜻하지 않음. WGS-only의 최적 메모리를 확정한 값도 아님.
- 클라우드 기본은 샘플당 고정16CPU임. HiFi Stage10은 이 예산 안에서 HAP을 순차 처리하며, 미지정 sort 버퍼만 RAM에 맞춰 최대8GiB로 확대함. 부족하면4GiB·기존값으로 낮춤.32CPU를 자동 요청하지 않음.
- 여러 샘플을 처리하려면 Docker CPU·RAM과 JSON의 풀/동시수를 함께 올림. 예를 들어 CPU68풀에서 샘플16CPU×4개를 설정 가능함. 실제 노드/다른 작업과 공유하는 커널 한도를 우선함.
- 서로 다른 노드/독립 컨테이너의 PID·CPU lease를 같은 state 폴더에 공유하지 않음. 노드별 state와 고유 run_id를 사용하고 동일 샘플 output 중복 실행을 피함.
- 이 구성의 cgroup 탐지 기준은 Linux cgroup v2임. v1 호스트에서 같은 한도 탐지가 된다고 검증한 것은 아님.

## 4. 중단과 재개

Tini가 PID1에서 종료 신호를 전달하며 entrypoint는 기존 run_pipeline.sh를 exec함. 별도 큐 제어기는 추가하지 않음.

```bash
# 실행한 컨테이너를 다른 터미널에서 중지함.
docker stop --timeout 120 mfm-hg00308
```

컨테이너를 제거해도 bind mount한 data/state/output은 남음. 원래 입력·작업파일·상태가 있어야 단계 재개가 가능하며 Docker가 삭제된 중간자료를 복원해주지는 않음.

위 단일샘플 실행을 Stage01부터 복구하는 예시는 다음과 같음.

```bash
docker run "${MFM_DOCKER_ARGS[@]}" mfm-preprocess:261124-cpu \
  resume --effective-config /data/state/cloud_azure_poc_HG00308/effective_config.json \
  --selection /data/state/cloud_azure_poc_HG00308/selected_samples.tsv \
  --sample HG00308 --from-stage 01
# 복구 요청 뒤 원래 명령으로 큐를 다시 시작함.
docker run "${MFM_DOCKER_ARGS[@]}" mfm-preprocess:261124-cpu \
  run --target azure_poc --sample HG00308 --config /config/job.json
```

전체 target 큐와 단일샘플 작업의 상태폴더는 다름. 전체 큐를 실행했다면 기존 문서의 `resume --target ... --sample ... --config ...`를 사용함.

## 5. 이미지 검사와 전달

```bash
# 실제 실행 이미지의 해시/환경 정보를 기록함.
docker image inspect --format '{{.Id}}' mfm-preprocess:261124-cpu
docker run --rm --entrypoint /bin/cat mfm-preprocess:261124-cpu \
  /opt/mfm-image/conda-explicit.txt

# 이미지 안에서 native·소형 WGS/RNA 검사10개를 다시 실행함.
docker run --rm --cpus 2 --memory 8g \
  --entrypoint /opt/conda/envs/jj_mfm/bin/python mfm-preprocess:261124-cpu \
  /opt/mfm/docker/smoke.py

# 레지스트리 없이 다른 서버로 전달할 때 사용함.
docker save mfm-preprocess:261124-cpu | gzip > mfm-preprocess_261124-cpu.tar.gz
# 대상 서버에서 실행함.
gunzip -c mfm-preprocess_261124-cpu.tar.gz | docker load
```

- 이미지의 /opt/mfm-image에는 실제 Conda explicit 목록·pip freeze·OS 패키지 목록이 남음. 빌드에서 검증한 native 바이너리는 첫 실행 때 쓰기 가능한 캐시에 잠금/원자적 복사로 재사용하며 기존 캐시를 덮어쓰지 않음.
- 소형 검사와 실제 donor 전과정의 처리시간·RAM·I/O는 구분함. 실제 클라우드 노드에서 한 샘플을 완료해 결과와 자원 사용량을 대조하는 PoC가 필요함.
- 현재 구현한 이미지는 CPU WGS/RNA 범위임. Clair3 self-call/CuPy/Parabricks GPU runtime은 포함하지 않음. 공개 VCF 기반 기본3212명 구성과 구분하며 GPU는 별도 이미지·호스트 runtime·경로 연동을 검증해야 함.
- 작성 서버에는 Docker 실행기가 없고 rootless namespace도 제한되어 로컬 image build/run은 수행하지 않음. entrypoint·기존 파이프라인 연결·native13종 및 소형검사는 별도 설치된 jj_mfm으로 검증함.
- [GitHub 자동 검사](../.github/workflows/docker.yml)는 Ubuntu24.04 실행기에서 실제 image build 후 컨테이너의 help/doctor/dry-run·소형 WGS/RNA 검사를 수행함. 해당 커밋의 PR Checks에서 성공 여부를 확인함. 검사 중 실제 donor 원본 다운로드·대용량 전처리·레지스트리 이미지 게시를 수행하지 않음.
- 자동 검사의 실행 CPU는2이며 설정의 샘플/풀 CPU도2로 맞춤. 제공 실사용 JSON의16CPU·384GiB 예제를 변경하지 않음. 작은 검사를 통과한 것만으로 해당 작은 자원에서 실제 donor 전처리가 가능하다고 판단하지 않음.

## 공식 문서

- [빌드와 환경 고정](https://docs.docker.com/build/building/best-practices/)
- [호스트 데이터 연결](https://docs.docker.com/engine/storage/bind-mounts/)
- [CPU·메모리 제한](https://docs.docker.com/engine/containers/resource_constraints/)
