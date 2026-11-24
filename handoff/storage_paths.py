"""클라우드 복사본의 명시 JSON 경로를 보존하고 상위 저장소 자동 연결을 막음."""
from __future__ import annotations
from pathlib import Path


def storage_defaults(config: dict, package: Path) -> dict:
    """상대 경로는 기존 load_config가 설정 파일 위치 기준으로 해석함."""
    return dict(config)
