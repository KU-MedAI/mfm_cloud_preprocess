# 완료 결과의 NAS 자동·수동 이관

NAS 또는 NAS가 마운트된 받는 서버에서 실행함. 전처리 큐는 기존 `run_pipeline.sh run`으로 계속 실행하며, 같은 수신 도구가 완료 결과를 자동 조회하거나 수동 배치로 이관함. 송신 서버에 별도 감시기·전처리 제어기는 추가하지 않음.

## 저장 위치

```text
/NFS_DATA3/mfm_data/preprocessed_data/
├── final/                         기존 서열/readset 데이터
│   ├── long_read_sequence/
│   └── short_read_sequence/
└── final_pileup/                  새 pileup 결과
    ├── short_read_pileup/
    │   ├── wgs_only/<project>/<sample>/
    │   └── pair/<project>/<sample>/
    └── long_read_pileup/
        ├── wgs_only/<project>/<sample>/
        └── pair/<project>/<sample>/
```

원본 서버에서 계산 중인 출력 경로는 유지함. NAS 쪽 `final_pileup`은 수신 명령이 생성함. 기존 `final` 자료와 합치거나 기존 폴더를 변경하지 않음.

manifest의 `read_type`으로 short/long을, `data_group=wgs` 또는 `wgs_rna_pair`로 wgs_only/pair를 구분함. 프로젝트는 `project_id`, 마지막 폴더는 `sample_id`임. 실행명·날짜는 NAS 데이터 경로에 넣지 않음. pair 한 샘플의 DNA와 RNA는 같은 샘플 폴더 아래 `compact/`·`rna_personal/`에 함께 보관함.

예: `final_pileup/short_read_pileup/pair/1KGP_ILLUMINA/HG00186/`임. 기존 준비 배치가 있으면 저장된 목적지로 재개하며 이 변경이 이미 전송된 폴더를 자동 이동하지 않음. 같은 프로젝트/샘플 목적지에 다른 자료가 있으면 덮어쓰지 않고 검증 오류로 중단함.

## 받는 서버에서 자동 실행

받는 서버에는 Python3·rsync·SSH client가 필요함. `--source-python`은 SSH 원본 서버에서 직접 실행 가능한 conda jj_mfm의 Python 경로임. Docker 전용 VM에서 host SSH와 컨테이너 실행을 연결하는 배포 구성은 별도임. 이미지에는 원본 명령과 rsync/SSH client를 포함함.

먼저 이 저장소의 `tools/pull_completed.py` 한 파일을 받는 서버로 복사함. 아래 예제의 SSH 주소·포트·키와 코드/상태/Python 경로를 실제 원본에 맞춤.

```bash
python3 pull_completed.py \
  --source USER@SOURCE_HOST --port SOURCE_SSH_PORT \
  --identity "$HOME/.ssh/korea_key" \
  --source-code /path/to/process_code \
  --source-config /path/to/state/effective_config.json \
  --source-python /path/to/conda/envs/jj_mfm/bin/python \
  --destination /NFS_DATA3/mfm_data/preprocessed_data/final_pileup \
  --watch --batch-size 1 --poll-interval 300 --remove-source
```

- 한 번 실행한 뒤 계속 조회함. 새 완료본이 없으면 기본300초(5분) 뒤 다시 확인하며 완료1명만 있어도 수신함. 수신이 끝나면 대기 없이 다음 완료본을 확인함.
- `--batch-size 1`은 새 완료본을 한 명씩 검증·서버 회수하기 위한 값임. 사용자가 샘플 이름을 지정할 필요는 없음. 기존 진행 중 배치가 있으면 그 배치부터 재개함.
- 이미 수신·검증하여 송신 서버에 이관 기록이 있는 샘플은 다음 선별에서 제외함. NAS 최종 폴더만 있고 송신 확인이 끝나지 않은 경우에는 내용을 대조하고 확인 단계를 재개하며 payload를 다시 복사하지 않음.
- 동일 NAS 목적지의 수신 잠금은 자동 대기 중에도 유지하여 두 수신기가 동시에 실행되지 않게 함. SSH/rsync 명령 실패는 다음 주기에 재시도하고 부분 전송을 이어받음. 수신측 checksum·파일 목록·기존 최종 경로 내용이 다르면 서버 사본을 보존하고 중단함.
- 자동 모드는 비대화 SSH를 사용함. 기존 키 또는 SSH agent로 접속할 수 있어야 함. 로그 시각은 수신 서버의 로컬 시간 기준이며 자동 모드에서는 전송 통계와 샘플 완료를 남김.
- SSH 터미널을 닫아도 계속 받으려면 tmux에서 실행함. 전체 선택 범위의 이관이 끝난 뒤에도 새 완료본 조회는 계속되며 중단은 Ctrl-C로 수행함. 현재 선택 범위는 지정한 `--source-config`의 샘플 목록임.

```bash
tmux new-session -s mfm_nas_pull
# 위 python3 명령을 실행한 뒤 Ctrl-B, D로 분리함.
```

## 수동10명 실행

자동 실행 명령에서 `--watch --poll-interval 300`을 빼고 `--batch-size 10`을 지정함. 완료10명 미만이면 현재 대기 인원만 출력하고 전송/자료 정리를 시작하지 않음. 자동 모드와 특정 배치 고정 `--batch ID`는 함께 사용할 수 없음.

## 공통 검증·회수 계약

- 완료/QC/정리 후 PASS·현재 출력 정책 일치·worker 잠금 해제된 샘플만 선택함. 이미 이관한 샘플은 다음 배치에서 제외함.
- 자동 실행 예제는1명씩, 수동 실행은10명을 한 배치로 받음. 배치 내부도 한 명씩 순서대로 받음. 현재 첫3명 기준 모델입력은1명 약76GB·10명 약760GB이며 부속 기록 용량은 별도임.
- rsync의 `--fsync`로 쓰기 완료를 요청하고 각 샘플의 모든 파일에 대해 SHA256을 대조함. `.incoming/<batch>/<sample>`에서 검증을 마친 뒤 최종 위치에 게시하고, 원본 서버에 수신 검증 결과를 전달함.
- `--remove-source`는 NAS 검증이 성공한 샘플에 한해 서버의 `compact/`·`rna_personal/`을 회수함. source 내용을 다시 대조하고 읽는 프로세스가 없는지 확인함. QC·메타데이터·로그·이관 영수증과 누적 완료 용량은 서버에 남음.
- 옵션을 빼면 NAS 검증·이관 기록까지 수행하고 서버 payload는 유지함. 추후 그 배치를 정리하려면 출력된 ID를 `--batch ID --remove-source`로 지정함.
- 중단·오류 시 같은 명령으로 재개함. 진행 중 배치 ID가 있으면 그 배치를 먼저 이어감. 이미 검증·게시한 NAS 샘플은 재복사하지 않고 내용 검증 후 진행함. 부분 전송의 같은 크기/mtime 내용 오류도 checksum 비교로 복구함.
- 기존 NAS 최종 경로에 다른 내용이 있으면 덮어쓰지 않고 중단함. 복사/검증 실패 샘플의 source payload는 보존함.

## 마지막 잔여 인원

자동 모드는 잔여 인원도 그대로 받음. 수동10명 모드에서는188명이10명씩18회와 마지막8명이며, 마지막 배치에 `--allow-partial`을 추가함.

## 완료·검증 기록

- 원본 서버: `<state_root>/nas_exports/<batch>/batch.json`과 샘플 `metadata/nas_archive.json`에 저장함.
- NAS: `final_pileup/_transfer_records/<batch>/`에 파일 목록과 검증 결과를 저장함.
- 이관 후에도 큐·재시작·worker는 완료 상태를 유지함. 이미 실행 중인 구 큐도 읽는 QC/완료 metadata를 보존하므로 다시 처리하지 않음.
- 카탈로그는 NAS 경로와 `processing_status=archived`를 기록함. source 호스트에서 NAS 파일 존재를 확인했다고 표시하지 않으며 `qc_all_paths_exist=False`, 로컬 `training_eligible=False`임. 이관 전 적격성은 `archive_training_eligible`로 보존함. 이관 자체가 생물학적 QC나 배열 값을 바꾸지 않음.
- NAS 전체 파일 검증 결과를 인증된 SSH로 원본에 전달하는 계약임. 도구 밖에서 샘플 폴더를 먼저 삭제하거나 영수증을 직접 작성하지 않음.
