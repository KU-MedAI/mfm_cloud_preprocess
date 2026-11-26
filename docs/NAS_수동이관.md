# 완료10명씩 NAS로 수동 이관

NAS 또는 NAS가 마운트된 받는 서버에서 실행함. 전처리 큐는 기존 `run_pipeline.sh run`으로 계속 실행하며, 이 도구는 완료 결과의 전송과 공간 회수만 수행함. 주기 실행기는 추가하지 않음.

## 저장 위치

```text
/NFS_DATA3/mfm_data/preprocessed_data/
├── final/                         기존 서열/readset 데이터
│   ├── long_read_sequence/
│   └── short_read_sequence/
└── final_pileup/                  새 pileup 결과
    └── <원본 output_root 폴더명>/
        └── wgs_rna_pair/<project>/<sample>/
```

원본 서버에서 계산 중인 출력 경로는 유지함. NAS 쪽 `final_pileup`은 수신 명령이 생성함. 기존 `final` 자료와 합치거나 기존 폴더를 변경하지 않음.

## 받는 서버에서 실행

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
  --batch-size 10 --remove-source
```

- 완료10명 미만이면 현재 대기 인원만 출력하고 전송/자료 정리를 시작하지 않음.
- 완료/QC/정리 후 PASS·현재 출력 정책 일치·worker 잠금 해제된 샘플만 선택함. 이미 이관한 샘플은 다음 배치에서 제외함.
- 10명은 한 배치이며, NAS 메모리·임시 공간을 줄이기 위해 한 명씩 순서대로 받음. 첫3명 기준 10명 모델입력은 약760GB이며 부속 기록 용량은 별도임.
- rsync의 `--fsync`로 쓰기 완료를 요청하고 각 샘플의 모든 파일에 대해 SHA256을 대조함. `.incoming/<batch>/<sample>`에서 검증을 마친 뒤 최종 위치에 게시하고, 원본 서버에 수신 검증 결과를 전달함.
- `--remove-source`는 NAS 검증이 성공한 샘플에 한해 서버의 `compact/`·`rna_personal/`을 회수함. source 내용을 다시 대조하고 읽는 프로세스가 없는지 확인함. QC·메타데이터·로그·이관 영수증과 누적 완료 용량은 서버에 남음.
- 옵션을 빼면 NAS 검증·이관 기록까지 수행하고 서버 payload는 유지함. 추후 그 배치를 정리하려면 출력된 ID를 `--batch ID --remove-source`로 지정함.
- 중단·오류 시 같은 명령으로 재개함. 진행 중 배치 ID가 있으면 그 배치를 먼저 이어감. 이미 검증·게시한 NAS 샘플은 재복사하지 않고 내용 검증 후 진행함. 부분 전송의 같은 크기/mtime 내용 오류도 checksum 비교로 복구함.
- 기존 NAS 최종 경로에 다른 내용이 있으면 덮어쓰지 않고 중단함. 복사/검증 실패 샘플의 source payload는 보존함.

## 마지막 잔여 인원

188명은10명씩18회, 마지막8명임. 마지막 배치만 위 명령에 `--allow-partial`을 추가함. 10명 미만도 명시적으로 수신함.

## 완료·검증 기록

- 원본 서버: `<state_root>/nas_exports/<batch>/batch.json`과 샘플 `metadata/nas_archive.json`에 저장함.
- NAS: `final_pileup/_transfer_records/<batch>/`에 파일 목록과 검증 결과를 저장함.
- 이관 후에도 큐·재시작·worker는 완료 상태를 유지함. 이미 실행 중인 구 큐도 읽는 QC/완료 metadata를 보존하므로 다시 처리하지 않음.
- 카탈로그는 NAS 경로와 `processing_status=archived`를 기록함. source 호스트에서 NAS 파일 존재를 확인했다고 표시하지 않으며 `qc_all_paths_exist=False`, 로컬 `training_eligible=False`임. 이관 전 적격성은 `archive_training_eligible`로 보존함. 이관 자체가 생물학적 QC나 배열 값을 바꾸지 않음.
- NAS 전체 파일 검증 결과를 인증된 SSH로 원본에 전달하는 계약임. 도구 밖에서 샘플 폴더를 먼저 삭제하거나 영수증을 직접 작성하지 않음.
