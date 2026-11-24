#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
TMP=$(mktemp -d /tmp/mfm-stage04-entry.XXXXXX)
trap '[[ "$TMP" == /tmp/mfm-stage04-entry.* ]] && rm -rf "$TMP"' EXIT

python - "$TMP" <<'PY'
import csv
import random
import sys
from pathlib import Path
import pysam
root = Path(sys.argv[1])
rng = random.Random(91)
sequence = "".join(rng.choice("ACGT") for _ in range(6000))
ref = root / "ref.fa"
ref.write_text(">chr1\n" + sequence + "\n")
pysam.faidx(str(ref))
ubam = root / "reads.ubam"
with pysam.AlignmentFile(ubam, "wb", header={"HD": {"VN": "1.6"}}) as out:
    for index, quality in enumerate(("I", "H")):
        read = pysam.AlignedSegment()
        read.query_name = f"read{index}"
        read.query_sequence = sequence[500 + index * 2000:1500 + index * 2000]
        read.query_qualities = pysam.qualitystring_to_array(quality * 1000)
        read.flag = 4
        out.write(read)
file_columns = ["sample_id", "run_id", "file_role", "format", "alignment_state", "path",
                "source_url", "size_bytes", "md5", "source_reference_fasta", "read_group_id"]
with (root / "files.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=file_columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerow({"sample_id": "S1", "run_id": "run1", "file_role": "alignment",
                     "format": "ubam", "alignment_state": "unaligned", "path": str(ubam),
                     "source_url": "-", "size_bytes": ubam.stat().st_size, "md5": "-",
                     "source_reference_fasta": "-", "read_group_id": "run1"})
sample_columns = ["sample_id", "platform", "read_type", "dataset_role", "project_id", "cohort",
                  "alignment_path", "raw_root", "rna_path", "final_root", "work_root",
                  "common_vcf_dir", "persample_vcf_dir", "rna_bigwigs", "vcf_source",
                  "public_vcf_path", "public_vcf_sample_id", "vcf_reference_fasta", "vcf_callset",
                  "phase_source", "read_files_manifest", "prepared_alignment_path",
                  "source_decode_reference_fasta", "analysis_reference_fasta"]
row = {key: "-" for key in sample_columns}
row.update({"sample_id": "S1", "platform": "hifi", "read_type": "long_read",
            "dataset_role": "pre_training", "project_id": "test_project", "cohort": "test_cohort",
            "work_root": str(root / "work"), "vcf_source": "clair3", "vcf_callset": "Clair3",
            "phase_source": "provided", "read_files_manifest": "files.tsv",
            "source_decode_reference_fasta": str(ref), "analysis_reference_fasta": str(ref)})
with (root / "samples.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=sample_columns, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerow(row)
PY

export PIPELINE_MANIFEST="$TMP/samples.tsv"
export DATASET_ROOT="$TMP/data"
export RUN_ID=stage04_test
export REF_DECODE="$TMP/ref.fa"
export REF_SEQUENCE="$TMP/ref.fa"
export AUTOSOMES=chr1
export CHROMS=chr1
export STAGE04_ALIGN_THREADS=2
export STAGE04_SORT_THREADS=1
export STAGE04_SORT_MEM=64M
export STAGE04_EXTRACT_THREADS=1
export STAGE04_MIN_QUAL_VALUES=2
export PREFLIGHT_MIN_FREE_GIB=0
export CLEAN_INTERMEDIATES_ON_SUCCESS=0
export ENABLE_CLAIR3_SELF_CALL=1

bash "$ROOT/01_prepare_analysis_bam.sh" S1
BAM="$TMP/work/prepared_input/S1.analysis.bam"
QC="$TMP/work/prepared_input/S1.stage04.json"
[[ -s "$BAM" && -s "$BAM.bai" && -s "$QC" ]]
bash "$ROOT/02_preflight.sh" S1
REPORT="$TMP/data/_run_status/stage04_test/wgs/test_project/S1/preflight/S1.preflight.tsv"
[[ -s "$REPORT" ]]
awk -F'\t' '$1=="stage04_reference_provenance" && $2=="PASS"{p=1}
             $1=="cram_decode_m5_chr1" && $2=="PASS"{m=1}
             END{exit(p&&m ? 0 : 1)}' "$REPORT"
echo "[PASS] Stage 04 entrypoint → Stage 05 provenance contract"

if WGS_REALIGN_FORMAT=bam bash -c "source '$ROOT/config_v2.sh'"; then
  echo '[FAIL] WGS_REALIGN_FORMAT=bam must be rejected before Stage 70' >&2
  exit 1
fi
