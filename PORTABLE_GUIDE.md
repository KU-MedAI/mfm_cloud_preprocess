# 클라우드용 외부 서버 안내

현재 복사본의 환경·자원·실행법은 [클라우드 실행](docs/클라우드_실행.md)과 [README](README.md)를 먼저 확인함. 아래 original의 고정68CPU 예시는 클라우드 기본값이 아니며 실제 할당과 JSON을 우선함.

# original 외부 서버 실행 안내

기준일: **2026-11-23**. 최신 모드·자원·정리 정책은 [GUIDE.md](GUIDE.md), 단계별 동작은 [PIPELINE.md](PIPELINE.md), 결과 구조는 [실행·출력 구조](docs/실행과_출력구조.md)를 기준으로 함.

## 1. 환경과 데이터 경로 준비

```bash
cd /path/to/cloud_전처리용
# Conda/Miniforge와 서버 기본 도구가 준비된 상태에서 환경 설치·보강함.
bash ./00_environment.sh
conda activate jj_mfm
bash ./00_environment.sh verify-core
# RNA를 처리하는 경우 확인함.
bash ./00_environment.sh verify-rna
cp handoff_config.example.json my_config.json
```

- `environment.yml` + `requirements.txt`로 `jj_mfm`, `clair3_environment.yml`로 선택 self-calling용 `clair3_v2`를 설치함. verify만 수행할 때는 인수 없는00번을 실행하지 않음.
- Bash/flock/GNU 도구·필요 시 Slurm은 서버에서 제공함. native pileup은 compiler/현재 conda HTSlib를 사용함. 실제 설치 버전은 `TOOL_VERSIONS.tsv`로 확인함.
- `my_config.json`의 `data_root`, `output_root`, `state_root`, `bundled_assets_root`, `run_id`를 새 서버에 맞춤. data/output/state는 서로 다른 경로임. 상대 경로는 JSON 파일 위치 기준임.
- 서버 전용 `handoff_config.json`·기존 `effective_config.json`을 경로 확인 없이 복사해 운영 실행하지 않음. 별도 실행은 output/state/run_id를 분리해 기존 완료본을 보존함.

## 2. payload와 자원 설정

- WGS `full`은 개인 HAP readset, `sequence_only`는 HAP 서열, `pileup_feature`는 HAP 서열 + reference pileup tensor임. JSON `pipeline_env`의 `MFM_SHORT_READ_OUTPUT_MODE`·`MFM_LONG_READ_OUTPUT_MODE`로 선택함. 제공 예제의 기본 출력은 `full`임.
- 기본 pileup WGS-only/coverage_only는 개인 HAP 재정렬을 하지 않음. `MFM_PILEUP_BENCHMARK_REALIGN=0`을 사용함. RNA `full`의 ASE/CNV용 WGS 근거 처리는 별개로 필요함.
- RNA `full`은 signal·ASE·CNV, `coverage_only`는 signal·QC임. RNA 보류는 `MFM_SEQUENCE_DEFER_RNA=1`이며 coverage_only와 동시에 사용하지 않음.
- 권장 자원 예시는 fixed_sample/샘플16CPU/동시4/풀68/RAM1.65TB·예약1.485TB임. 새 서버의 실제 CPU/RAM/Slurm/cgroup 한도에 맞춰 낮춰야 함. 고정 샘플 예산×동시수가 풀을 넘지 않게 설정함.
- 단독68CPU 시험은 sample_cpus68과local_jobs1을 함께 지정함. local_jobs1만 바꾸면 샘플16CPU가 그대로임. 자원 배정과 실제 사용량은 별도 측정함.
- 데이터 정리 전까지 필요한 공간은 최종 모델 용량보다 큼. 입력·분석 BAM·haplotag·정렬 임시·RNA STAR index 등을 포함해 계산함. `min_output_free_bytes`로 신규 입장 공간 조건을 지정 가능함.

## 3. 공식 내장 큐 실행

새 클라우드 분담은 `--target azure_poc|azure_main|aws_main|server`를 사용함. 예: `bash run_pipeline.sh run --target azure_poc --sample HG00308`임. target 명단/설정을 자동 선택하며 기존 lab/batch/selection을 함께 넣지 않음. 전체 명단은 CPU/RAM 예산 안에서 순차·병렬 처리하고 기존220명/10TB 조기중지를 적용하지 않음.

아래는 기존 lab/batch 방식의 호환 예시임.

```bash
# 한 기관·한 배치부터 확인함. 아래 selection은 실제 제공된 목록을 선택함.
bash run_pipeline.sh dry-run   --selection manifests/리드별/medai/001/long_read.tsv   --config my_config.json --executor local

# 최종 경로·payload·자원·다운로드 목록을 확인한 뒤 실제 실행함.
bash run_pipeline.sh run   --selection manifests/리드별/medai/001/long_read.tsv   --config my_config.json --executor local
```

- 내장 local 큐가 선행 다운로드·준비 완료 샘플 즉시 투입·자원 예약·검증·정리를 이어감. 별도 연속 제어 Python이 필수인 구조가 아님.
- `--executor slurm`은 Slurm config로 제출하며 실제 할당 범위를 우선함. 서로 다른 노드의 CPU 풀을 공유 잠금 하나로 묶지 않음.
- 기존 lab/batch 선택의 한 실행은 한 기관·한 배치임. 모든3,212명을 한 선택 TSV로 넘기지 않음. 새 분담 전체는 target으로 지정함. 여러 배치는 `run_read_group.sh` 또는 검토한 배치별 실행 목록으로 진행함.
- 순수 직접 branch 실행은 입력을 이미 준비한 경우에 사용함. `PIPELINE_MANIFEST`·INPUT/OUTPUT/자원 환경을 명시해야 하며 내장 다운로드 큐와 같은 기능으로 간주하지 않음.

## 4. reference와 입력 검증

Reference FASTA/FAI·VCF/index·GENCODE는 도구 설치와 별개 데이터임. 실제 모드·프로젝트에 필요한 파일은 dry-run의 입력 목록과 manifest로 결정함.

- Illumina CRAM은 정확히 일치하는 decode reference가 필요함. 파일 이름만 GRCh38이라고 바꾸지 않음.
- HGSVC/HPRC 등 프로젝트별 analysis/VCF reference 차이를 유지함. Stage02는 chr1–22 SQ/LN/M5(존재 시)와 reference/VCF donor·REF 계약을 확인함.
- RNA는 GENCODE v47 GTF와 GTEx 추가 contig reference·FAI가 필요함. index fingerprint를 확인하며 다른 HAP/프로필의 STAR index를 임의 재사용하지 않음.
- Blacklist는 기본 정책에서 창 삭제용이 아닌 주석용임. 선택 목록과 설정에 따라 요구되는 입력을 확인함. reference/VCF 불일치를 검증 우회로 통과시키지 않음.
- 한 donor의 여러 run/uBAM/FASTQ는 manifest에 연결된 모든 입력을 사용함. 한 파일만 임의 선택하거나 별개 donor로 집계하지 않음.

## 5. 최종 출력과 정리

- 결과는 `<배치 결과 root>/<wgs|wgs_rna_pair>/<project>/<sample>/` 아래이며 JSON `output_layout`에 따라 배치 root가 결정됨.
- 최종 모델은 `compact/`와 해당 RNA Zarr/TSV임. 주석·압축 chain·QC·metadata·로그를 함께 전달함. 원본/중간/계측 자료의 용량과 구분함.
- 최종QC→삭제전 guard→중간 정리→정리후 검증 순서를 따름. 기존 원본을 보존해야 하는 경우 `remove_owned_raw_inputs=false`뿐 아니라 `MFM_MODEL_ONLY_CLEANUP=0` 등 [GUIDE의 조합](GUIDE.md#7-검증재개정리)을 사용함.
- 중간 정리 완료 뒤에는 개인 FASTA/BAM/STAR index가 없을 수 있음. RNA 추가·재개는 같은 reference/VCF/window 계약을 재생성·검증한 뒤 수행함.
- 실패는 샘플별 상태와 로그에 남음. `resume --sample ID --from-stage NN`은 검토된 WGS 단계의 내장 재개 요청임. 전체 큐가 끝난 상태이면 같은 config·selection의 `run`도 다시 실행함.

## 6. 코드 전달·기존 shard 자료

전달 bundle을 새 경로에 만들 때만 아래를 실행함. 기존 경로는 덮어쓰지 않음.

```bash
bash scripts/make_handoff_bundle.sh /new/path/cloud_bundle
cd /new/path/cloud_bundle
python3 tools/check_cloud_package.py
```

- `share_package/`는 원본 출처·동결 registry 근거임. 파일 수는 인원수와 다르며 공통 reference/panel/annotation과 샘플 전용 read/VCF/index를 구분함.
- 2026-09-28의37-shard·10TiB·ETA·단계별4→4→2 예시는 당시 readset 처리 가정임. 현재 기본 큐/모드의 동시수·용량·기간으로 그대로 사용하지 않음. 기존 방식 재현은 [SHARDED_HANDOFF.md](SHARDED_HANDOFF.md)를 참조함.
- 전달 manifest의 library provenance는 [LIBRARY_PROVENANCE.md](LIBRARY_PROVENANCE.md)를 따름. `pcr_free_inferred`를 `confirmed`로 자동 변경하지 않음.
- 기존37-shard snapshot을 사용하는 경우 아래 확인을 수행함. 해당 snapshot의 정상 기준은 `shards=37 samples=3212 mode=check`이며 신규 shard에도 적용되는 고정 숫자는 아님.

```bash
python scripts/sync_shard_provenance.py --check-only
```
