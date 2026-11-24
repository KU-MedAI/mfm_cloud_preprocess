"""RNA 재사용을 위한 입력·설정·완료 출력 계약. 대형 입력의 반복 전체 스캔을 피한다."""
from __future__ import annotations

import hashlib
import json
import shutil
import time
from pathlib import Path
from typing import Any

CACHE_SCHEMA = "mfm_rna_cache_v1"


def _stat(path: Path) -> dict[str, int]:
    value = path.stat()
    return {"size_bytes": value.st_size, "mtime_ns": value.st_mtime_ns,
            "ctime_ns": value.st_ctime_ns, "inode": value.st_ino, "device": value.st_dev}


class FingerprintChanged(ValueError):
    """지문을 읽는 동안 파일 버전이 바뀜."""


def file_fingerprint(path: str | Path, *, content: bool = False) -> dict[str, Any]:
    for attempt in range(3):
        try:
            return _file_fingerprint_once(path, content=content)
        except FingerprintChanged:
            if attempt == 2:
                raise
            time.sleep(.2)


def _file_fingerprint_once(path: str | Path, *, content: bool = False) -> dict[str, Any]:
    """계약 파일은 전체 SHA256, 큰 서열 파일은 stat와 내용 표본으로 변경을 감지한다."""
    path = Path(path).resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"캐시 입력이 파일이 아님: {path}")
    before = _stat(path)
    result = {"type": "file", "path": str(path), **before}
    if content:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            while block := handle.read(8 << 20):
                digest.update(block)
        result["sha256"] = digest.hexdigest()
    else:
        # 같은 파일시스템 clock tick 안의 같은 크기 덮어쓰기도 표본 구간에서 확인한다.
        # 1MiB 이하는 전체를 읽고 큰 파일은 최대 8×64KiB를 읽는다.
        digest = hashlib.sha256()
        size = before["size_bytes"]
        width = 1 << 16
        with path.open("rb") as handle:
            if size <= 1 << 20:
                digest.update(handle.read())
            else:
                for offset in sorted({(size - width) * i // 7 for i in range(8)}):
                    handle.seek(offset)
                    digest.update(str(offset).encode() + b":" + handle.read(width))
        result["probe_sha256"] = digest.hexdigest()
        result["probe_policy"] = "full_le_1MiB_else_8x64KiB"
    if _stat(path) != before:
        raise FingerprintChanged(f"fingerprint 계산 중 파일 변경: {path}")
    return result


def optional_file(path: str | Path | None, *, content: bool = False) -> dict | None:
    return file_fingerprint(path, content=content) if path else None


def tool_fingerprint(command: str) -> dict:
    path = shutil.which(command)
    if path is None:
        raise ValueError(f"실행 도구 없음: {command}")
    return file_fingerprint(path)


def fingerprint(payload: Any) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=False).encode()).hexdigest()


def files_unchanged(value: Any) -> bool:
    """계약에 들어 있는 파일을 재확인한다. 생성 도중 입력이 바뀌면 완료하지 않는다."""
    try:
        if isinstance(value, dict):
            if value.get("type") == "file":
                return file_fingerprint(value["path"], content="sha256" in value) == value
            return all(files_unchanged(item) for item in value.values())
        if isinstance(value, (list, tuple)):
            return all(files_unchanged(item) for item in value)
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def cache_fields(inputs: dict, outputs: list[Path]) -> dict:
    if not files_unchanged(inputs):
        raise ValueError("RNA 산출물 생성 중 입력이 변경됨; 완료 cache 기록 중단")
    return {"cache_schema": CACHE_SCHEMA, "cache_fingerprint": fingerprint(inputs),
            "cache_inputs": inputs, "cache_outputs": [file_fingerprint(p) for p in outputs]}


def cache_valid(record: Any, inputs: dict | None = None) -> bool:
    if not isinstance(record, dict) or record.get("cache_schema") != CACHE_SCHEMA:
        return False
    stored = record.get("cache_inputs")
    outputs = record.get("cache_outputs")
    return (isinstance(stored, dict) and isinstance(outputs, list) and bool(outputs)
            and all(isinstance(item, dict) and item.get("type") == "file" for item in outputs)
            and record.get("cache_fingerprint") == fingerprint(stored)
            and (inputs is None or fingerprint(stored) == fingerprint(inputs))
            and files_unchanged(stored) and files_unchanged(outputs))


def read_cache(path: str | Path) -> dict:
    try:
        result = json.loads(Path(path).read_text())
        return result if isinstance(result, dict) else {}
    except (OSError, ValueError):
        return {}
