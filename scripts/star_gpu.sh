#!/usr/bin/env bash
set -euo pipefail
GPU_PACKAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="$GPU_PACKAGE_DIR/src:${PYTHONPATH:-}"
exec "${PYTHON_BIN:-python3}" -B -m integrate_wgs_rna.gpu_rna "$@"
