"""자체 생성 검증 기록을 디스크에 순차 저장하고 반복 재생한다."""
from pathlib import Path
import pickle


class ReplayRecords:
    def __init__(self, path: Path):
        self.path = path
        self.handle = path.open('wb', buffering=1 << 20)
        self.count = 0

    def append(self, record: dict):
        pickle.dump(record, self.handle, protocol=5)
        self.count += 1

    def __iter__(self):
        if self.handle is not None:
            self.handle.flush()
        with self.path.open('rb', buffering=1 << 20) as source:
            for _ in range(self.count):
                yield pickle.load(source)

    def close(self):
        if self.handle is not None:
            self.handle.close()
            self.handle = None
