#!/usr/bin/env bash
#SBATCH --job-name=mfm_wgs_v2
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-100%4
#SBATCH --cpus-per-task=64
#SBATCH --mem=512G
#SBATCH --time=96:00:00

# integrated v2 WGS SLURM array.
# 사용 예:
#   sbatch --array=1-20%4 --export=ALL,SAMPLE_LIST=/path/samples.txt slurm_wgs.sbatch
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config_v2.sh"
LIST=${SAMPLE_LIST:?SAMPLE_LIST를 지정하세요. 절대경로 권장}
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[[ -n "$S" ]] || { echo "empty sample at array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
export CLAIR3_THREADS=${SLURM_CPUS_PER_TASK:-32}
export WGS_REALIGN_THREADS=${SLURM_CPUS_PER_TASK:-32}
export COMPACT_JOBS=${COMPACT_JOBS:-${SLURM_CPUS_PER_TASK:-32}}
export DENOVO_THREADS=${DENOVO_THREADS:-${SLURM_CPUS_PER_TASK:-64}}
export ASSEMBLY_ANCHOR_THREADS=${ASSEMBLY_ANCHOR_THREADS:-${SLURM_CPUS_PER_TASK:-64}}
echo "[$(date '+%F %T')] WGS $S array=${SLURM_ARRAY_TASK_ID}"
exec "$SCRIPT_DIR/run_wgs_sample.sh" "$S"
