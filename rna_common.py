"""RNA 단계 스크립트(13, 15–17, 19) 공통 도우미 — 로그, 외부 도구 실행, 원자적 JSON 쓰기."""
from __future__ import annotations

import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def log(msg: str) -> None:
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", file=sys.stderr, flush=True)


def warn(msg: str) -> None:
    log(f"[WARN] {msg}")


def die(msg: str, code: int = 1) -> "None":
    log(f"[FATAL] {msg}")
    raise SystemExit(code)


def run(cmd: list[str], *, log_path: Path | None = None, what: str = "") -> None:
    """외부 명령 실행. 실패하면 로그 마지막 3줄과 함께 중단."""
    if os.environ.get('MFM_CPU_POOL_LEASE_FILE'):
        sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
        from integrate_wgs_rna.rna_memory import prepare_star_resources
        prepare_star_resources(cmd)
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("w") as fh:
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
    else:
        rc = subprocess.run(cmd).returncode
    if rc != 0:
        tail = ""
        if log_path is not None and log_path.is_file():
            tail = " ".join(log_path.read_text(errors="replace").splitlines()[-3:])
        die(f"{what or cmd[0]} 실패 (exit {rc}) {tail}")


def write_json(path: Path, payload: dict) -> None:
    """.part 에 쓴 뒤 원자적으로 교체한다 — 중간에 죽어도 반쪽 파일이 남지 않는다."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    os.replace(tmp, path)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def nonempty(*paths: Path) -> bool:
    return all(p.is_file() and p.stat().st_size > 0 for p in paths)
