#!/usr/bin/env bash
#SBATCH --job-name=mfm_rna
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-100%6            # STAR index 가 샘플당 26 GB 쓰므로 동시 실행수를 디스크로 제한할 것
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G                  # STAR index/정렬 메모리(2.8 Gb 게놈 기준 여유 포함)
#SBATCH --time=12:00:00
#
# RNA 파이프라인 slurm array. WGS(개인 서열)가 끝난 뒤 제출할 것.
#   sbatch --dependency=afterok:<WGS_JOBID> --export=ALL,SAMPLE_LIST=/path/list.txt slurm_rna.sbatch
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBMIT_DIR=${SLURM_SUBMIT_DIR:-$PWD}
LIST=${SAMPLE_LIST:?SAMPLE_LIST 를 지정할 것}
[[ "$LIST" == /* ]] || LIST="$SUBMIT_DIR/$LIST"
[[ -s "$LIST" ]] || { echo "sample list 없음: $LIST" >&2; exit 1; }
cd "$SCRIPT_DIR"
source "$SCRIPT_DIR/config_v2.sh"
S=$(awk -v target="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID 없음}" 'NF && $1 !~ /^#/ {n++; if (n == target) {print $1; exit}}' "$LIST")
[[ -n "$S" ]] || { echo "유효 sample 없음: index ${SLURM_ARRAY_TASK_ID} list=$LIST" >&2; exit 1; }
export STAR_INDEX_THREADS=${SLURM_CPUS_PER_TASK:-16}
export STAR_ALIGN_THREADS=${SLURM_CPUS_PER_TASK:-16}
echo "[$(date '+%F %T')] RNA $S (array ${SLURM_ARRAY_TASK_ID})"
exec "$SCRIPT_DIR/run_rna_sample.sh" "$S"
