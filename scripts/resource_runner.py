#!/usr/bin/env python3
"""명령 하나의 wall/CPU/max-RSS를 독립 JSON으로 기록한다."""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--address-space-bytes", type=int, default=0)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("-- 뒤에 실행할 명령이 필요합니다")
    if args.address_space_bytes < 0:
        parser.error('주소 공간 제한은 음수일 수 없습니다')
    if args.address_space_bytes:
        _, previous = resource.getrlimit(resource.RLIMIT_AS)
        target = args.address_space_bytes if previous == resource.RLIM_INFINITY else min(previous, args.address_space_bytes)
        resource.setrlimit(resource.RLIMIT_AS, (target, previous))

    boundary_info = None
    if os.environ.get('MFM_CPU_POOL_LEASE_FILE'):
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
        from integrate_wgs_rna.realign.boundary_upgrade import prepare_boundary
        command,boundary_info=prepare_boundary(command,args.output)
    started = time.time()
    monotonic_started = time.perf_counter()
    completed = subprocess.run(command, check=False)
    if boundary_info is not None and completed.returncode == 0:
        from integrate_wgs_rna.realign.boundary_upgrade import finish_boundary
        finish_boundary(boundary_info)
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    payload = {
        "schema": "mfm_resource_usage_v1",
        "command": command,
        "exit_code": completed.returncode,
        "started_at_unix": started,
        "elapsed_sec": round(time.perf_counter() - monotonic_started, 6),
        "user_cpu_sec": round(usage.ru_utime, 6),
        "system_cpu_sec": round(usage.ru_stime, 6),
        "max_rss_kib": int(usage.ru_maxrss),
        "address_space_limit_bytes": args.address_space_bytes,
    }
    if boundary_info is not None:payload['boundary_upgrade']=boundary_info
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.part.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, args.output)
    return completed.returncode


if __name__ == "__main__":
    sys.exit(main())
