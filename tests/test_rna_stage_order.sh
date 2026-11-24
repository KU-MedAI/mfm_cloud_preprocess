#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUNNER=$ROOT_DIR/run_rna_sample.sh
bash -n "$RUNNER"

python3 - "$RUNNER" <<'PY'
import re,sys
text=open(sys.argv[1]).read()
positions=[]
for stage in range(13,27):
    match=re.search(rf'(?:pipeline_(?:run_command|run_stage|skip_stage|run_part)|"\$RUNNER") {stage} ', text)
    assert match, f"stage {stage} orchestration missing"
    positions.append(match.start())
assert positions == sorted(positions), positions
assert '--preflight-json "$F/qc/rna_preflight.json"' in text
assert '--windows-tsv "$F/windows/windows.tsv"' in text
assert '--window-out "$TD/coverage_windows.zarr"' in text
assert '25_validate_sample_package.sh "$S" rna' in text
assert '26_cleanup_intermediates.sh rna' in text
assert '--hap-variants "$HV"' in text
PY

grep -Fq 'RNA 13–26, 검증 25, cleanup 26' "$ROOT_DIR/run_sample.sh"

echo "[PASS] RNA 13-26 numbered order/preflight/dual-output-coverage/validated-HAP/finalize contract"
