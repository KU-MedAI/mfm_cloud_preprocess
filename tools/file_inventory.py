"""검사 시점의 파일 stat을 한 번씩 읽고 링크 디렉터리는 내려가지 않음."""
import os
from pathlib import Path


def file_stats(root):
    root = Path(root)
    if not root.is_dir():return
    pending = [root]
    while pending:
        with os.scandir(pending.pop()) as entries:
            for entry in entries:
                if entry.is_dir(follow_symlinks=False):
                    pending.append(entry.path)
                elif entry.is_file():
                    # 기존 Path.is_file/stat처럼 파일 링크의 실제 대상도 검사함.
                    yield Path(entry.path), entry.stat()
