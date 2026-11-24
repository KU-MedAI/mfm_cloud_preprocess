#!/usr/bin/env python3
"""실패 단계의 원본 로그를 연결한다. 자동 추출 문구를 확정 원인으로 오인하지 않는다."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from stage_provenance import append_jsonl


def main() -> int:
    parser = argparse.ArgumentParser()
    for name in ("output", "sample", "branch", "stage", "name", "log"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--run-id", default="")
    parser.add_argument("--command-json", default="[]")
    parser.add_argument("--exit-code", type=int, required=True)
    args = parser.parse_args()
    path = Path(args.log)
    with path.open("rb") as handle:
        handle.seek(max(0, path.stat().st_size - 65536))
        tail = handle.read().decode("utf-8", errors="replace")
    signatures = [line for line in tail.splitlines()
                  if any(token in line.lower() for token in ("error", "fatal", "[fail]", "exception"))]
    payload = {
        "schema": "mfm_stage_failure_v1", "timestamp": datetime.now(timezone.utc).isoformat(),
        "sample_id": args.sample, "branch": args.branch, "stage": args.stage,
        "stage_name": args.name, "run_id": args.run_id, "exit_code": args.exit_code,
        "command_argv": json.loads(args.command_json), "log_path": str(path.resolve()),
        "log_bytes": path.stat().st_size, "error_signature": "\n".join(signatures[-12:]),
        "diagnosis_status": "unclassified", "log_tail": tail[-16384:],
    }
    append_jsonl(Path(args.output), payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
