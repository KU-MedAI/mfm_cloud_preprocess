#!/usr/bin/env python3
"""샘플별 stage 실행 provenance JSONL writer."""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def tool_version(path: str) -> str:
    if not path:
        return ""
    candidates = ([path, "--version"], [path, "version"])
    for command in candidates:
        try:
            result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, timeout=15)
        except (OSError, subprocess.SubprocessError):
            continue
        text = (result.stdout or "").strip().splitlines()
        if text:
            return text[0][:500]
    return ""


def append_jsonl(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.write(json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("event", choices=["start", "finish", "skip"])
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--branch", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--run-id", default="")
    parser.add_argument("--command-json", default="[]")
    parser.add_argument("--status", default="")
    parser.add_argument("--exit-code", type=int)
    parser.add_argument("--elapsed-seconds", type=float)
    parser.add_argument("--reason", default="")
    parser.add_argument("--log", default="")
    parser.add_argument("--qc", default="")
    parser.add_argument("--tool", action="append", default=[])
    args = parser.parse_args()

    try:
        command = json.loads(args.command_json)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"invalid --command-json: {exc}")
    tools: dict[str, dict[str, str]] = {}
    if args.event == "start":
        for spec in args.tool:
            name, _, path = spec.partition("=")
            if name and path:
                tools[name] = {"path": path, "version": tool_version(path)}

    payload = {
        "schema": "mfm_execution_provenance_v1",
        "event": args.event,
        "timestamp": datetime.now(timezone.utc).astimezone().isoformat(),
        "sample_id": args.sample,
        "branch": args.branch,
        "stage": args.stage,
        "stage_name": args.name,
        "run_id": args.run_id,
        "pid": os.getpid(),
        "command_argv": command,
    }
    if tools:
        payload["tools"] = tools
    if args.status:
        payload["status"] = args.status
    if args.exit_code is not None:
        payload["exit_code"] = args.exit_code
    if args.elapsed_seconds is not None:
        payload["elapsed_seconds"] = args.elapsed_seconds
    if args.reason:
        payload["reason"] = args.reason
    if args.log:
        payload["log_path"] = args.log
    if args.qc:
        payload["qc_path"] = args.qc
    append_jsonl(args.output, payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
