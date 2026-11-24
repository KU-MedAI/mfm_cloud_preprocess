#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SHARD_ROOT="$ROOT"
exec "$ROOT/pipeline/run_shard_queue.sh" "$@"
