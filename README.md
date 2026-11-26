# MFM 전처리 · Docker

Docker 엔진이 설치된 Linux x86_64 VM에서 Bash로 실행함. 기본 출력은 WGS pileup36채널과 RNA full(ASE/CNV 포함)임.

## 1. 이미지 빌드

```bash
git clone https://github.com/KU-MedAI/mfm_cloud_preprocess.git
cd mfm_cloud_preprocess
docker build --platform linux/amd64 -f docker/Dockerfile -t mfm-preprocess:261125-cpu .
```

## 2. 저장 경로·실행 설정

```bash
MFM_DOCKER_ROOT="$PWD/../mfm_docker_run"
mkdir -p "$MFM_DOCKER_ROOT"/data/{input,assets,output,state} "$MFM_DOCKER_ROOT"/{scratch,cache}
cp docker/config.pileup.example.json "$MFM_DOCKER_ROOT/config.json"
```

예제는 **샘플당 고정16CPU·384GiB·동시1명** 설정임. VM에 맞춰 RAM 예산과 아래 `--memory`를 함께 수정함. CPU 변경 시 JSON과 `--cpus`도 함께 맞춤.

```bash
MFM_DOCKER_ARGS=(
  --rm --user "$(id -u):$(id -g)"
  --cpus 16 --memory 384g --stop-timeout 120
  --mount "type=bind,src=$MFM_DOCKER_ROOT/data,dst=/data"
  --mount "type=bind,src=$MFM_DOCKER_ROOT/scratch,dst=/scratch"
  --mount "type=bind,src=$MFM_DOCKER_ROOT/cache,dst=/cache"
  --mount "type=bind,src=$MFM_DOCKER_ROOT/config.json,dst=/config/job.json,readonly"
)
```

## 3. 단일 샘플 전처리

```bash
# 사전 점검과 처리 계획 확인
docker run "${MFM_DOCKER_ARGS[@]}" mfm-preprocess:261125-cpu \
  doctor --config /config/job.json
docker run "${MFM_DOCKER_ARGS[@]}" mfm-preprocess:261125-cpu \
  dry-run --target azure_poc --sample HG00308 --config /config/job.json

# 실제 다운로드·전처리 실행
docker run "${MFM_DOCKER_ARGS[@]}" --name mfm-hg00308 mfm-preprocess:261125-cpu \
  run --target azure_poc --sample HG00308 --config /config/job.json
```

## 4. 대상 전체 전처리

```bash
docker run "${MFM_DOCKER_ARGS[@]}" mfm-preprocess:261125-cpu \
  run --target azure_poc --config /config/job.json
```

`--target`은 `azure_poc`, `azure_main`, `aws_main`, `server` 중 선택함. 선택한 대상 전체를 내장 큐로 처리함.

결과는 `$MFM_DOCKER_ROOT/data/output`, 실행 상태는 `$MFM_DOCKER_ROOT/data/state`에 저장함. 필요한 입력은 자동 수집하며 최종QC 뒤 중간파일·소유 원본을 정리함.

## 상세 설명

- [Docker 설정·중단/재개·이미지 전달](docs/Docker_실행.md)
- [환경·자원·출력 모드·전처리 정책](docs/전처리_상세.md)
- [완료 결과 NAS `final_pileup`으로 자동·수동 이관](docs/NAS_수동이관.md)

- [Pileup36채널·BAF·HP별 관측](docs/pileup_36채널.md)
