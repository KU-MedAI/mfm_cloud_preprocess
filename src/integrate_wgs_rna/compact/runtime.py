"""메모리 프로필이 조정하는 처리 전용 버퍼.

storage chunk/schema 상수와 분리한다. 이 값은 읽기·검증 batch 크기만 바꾸므로
동일 입력의 Zarr 배열 내용과 필터 정책을 바꾸지 않는다.
"""

from __future__ import annotations

import os


def positive_env_int(name: str, default: int) -> int:
    raw = os.environ.get(name, str(default))
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(f"{name} must be a positive integer: {raw}") from exc
    if value <= 0:
        raise ValueError(f"{name} must be a positive integer: {raw}")
    return value


def nonnegative_env_int(name: str, default: int) -> int:
    raw = os.environ.get(name, str(default))
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(f'{name}은 0 이상의 정수여야 합니다: {raw}') from exc
    if value < 0:
        raise ValueError(f'{name}은 0 이상의 정수여야 합니다: {raw}')
    return value


VALIDATION_RECORD_CHUNK = positive_env_int("MFM_VALIDATION_RECORD_CHUNK", 262_144)
VALIDATION_PACKED_BYTE_CHUNK = positive_env_int(
    "MFM_VALIDATION_PACKED_BYTE_CHUNK", 16 << 20
)
VALIDATION_CACHE_CHUNKS = nonnegative_env_int("MFM_VALIDATION_CACHE_CHUNKS", 4)

