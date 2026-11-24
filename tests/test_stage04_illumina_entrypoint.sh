#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
TMP=$(mktemp -d /tmp/mfm-stage04-illumina.XXXXXX)
trap '[[ "$TMP" == /tmp/mfm-stage04-illumina.* ]] && rm -rf "$TMP"' EXIT

python - "$TMP" <<'PY'
import csv
import gzip
import random
import sys
from pathlib import Path
import pysam

root = Path(sys.argv[1])
rng = random.Random(407)
sequence = "".join(rng.choice("ACGT") for _ in range(7000))
reference = root / "ref.fa"
reference.write_text(">chr1\n" + sequence + "\n")
pysam.faidx(str(reference))
r1, r2 = root / "R1.fastq.gz", root / "R2.fastq.gz"
with gzip.open(r1, "wt") as handle:
    handle.write(f"@pair1\n{sequence[1000:1150]}\n+\n{'I' * 150}\n")
with gzip.open(r2, "wt") as handle:
    handle.write(f"@pair1\n{sequence[1200:1350]}\n+\n{'H' * 150}\n")
file_columns = ["sample_id", "run_id", "file_role", "format", "alignment_state", "path",
                "source_url", "size_bytes", "md5", "source_reference_fasta", "read_group_id"]
with (root / "files.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=file_columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for role, path in (("read1", r1), ("read2", r2)):
        writer.writerow({"sample_id": "ILM1", "run_id": "runA", "file_role": role,
                         "format": "fastq_gz", "alignment_state": "unaligned", "path": str(path),
                         "source_url": "-", "size_bytes": path.stat().st_size, "md5": "-",
                         "source_reference_fasta": "-", "read_group_id": "rgA"})
sample_columns = ["sample_id", "platform", "read_type", "dataset_role", "project_id", "cohort",
                  "alignment_path", "raw_root", "rna_path", "final_root", "work_root",
                  "common_vcf_dir", "persample_vcf_dir", "rna_bigwigs", "vcf_source",
                  "public_vcf_path", "public_vcf_sample_id", "vcf_reference_fasta", "vcf_callset",
                  "phase_source", "read_files_manifest", "prepared_alignment_path",
                  "source_decode_reference_fasta", "analysis_reference_fasta"]
row = {key: "-" for key in sample_columns}
row.update({"sample_id": "ILM1", "platform": "illumina", "read_type": "short_read",
            "dataset_role": "pre_training", "project_id": "test_project", "cohort": "test_cohort",
            "work_root": str(root / "work"), "vcf_source": "illumina_panel", "vcf_callset": "1000G_panel",
            "phase_source": "provided", "read_files_manifest": "files.tsv",
            "source_decode_reference_fasta": str(reference), "analysis_reference_fasta": str(reference)})
with (root / "samples.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=sample_columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerow(row)
PY

export PIPELINE_MANIFEST="$TMP/samples.tsv"
export DATASET_ROOT="$TMP/data"
export RUN_ID=stage04_illumina_test
export REF_DECODE="$TMP/ref.fa"
export REF_SEQUENCE="$TMP/ref.fa"
export ILLUMINA_REF_DECODE="$TMP/ref.fa"
export AUTOSOMES=chr1
export CHROMS=chr1
export STAGE04_ALIGN_THREADS=2
export STAGE04_SORT_THREADS=1
export STAGE04_SORT_MEM=64M
export STAGE04_EXTRACT_THREADS=1
export STAGE04_MIN_QUAL_VALUES=2
export BWA_MEM2_INDEX_ROOT="$TMP/bwa_index_cache"
export PREFLIGHT_MIN_FREE_GIB=0

bash "$ROOT/01_prepare_analysis_bam.sh" ILM1
QC="$TMP/work/prepared_input/ILM1.stage04.json"
BAM="$TMP/work/prepared_input/ILM1.analysis.bam"
[[ -s "$QC" && -s "$BAM" && -s "$BAM.bai" ]]

python - "$QC" "$TMP/ref.fa" "$TMP/bwa_index_cache" <<'PY'
import json
import sys
from pathlib import Path

qc = json.loads(Path(sys.argv[1]).read_text())
reference = Path(sys.argv[2])
cache_root = Path(sys.argv[3]).resolve()
assert qc["status"] == "created_raw_illumina"
assert qc["alignment_policy"]["aligner"] == "bwa-mem2"
assert qc["preset"] == "mem"
assert qc["bwa_mem2_index"]["status"] == "created"
assert Path(qc["bwa_mem2_index"]["index_prefix"]).is_relative_to(cache_root)
assert qc["tools"]["bwa_mem2"].startswith("bwa-mem2 ")
assert qc["tools"]["picard"]
assert qc["duplicate_marking"]["action"] == "picard_marked"
assert qc["duplicate_marking"]["remove_duplicates"] is False
assert qc["duplicate_marking"]["before"]["records"] == 2
assert qc["duplicate_marking"]["after"]["records"] == 2
assert Path(qc["duplicate_marking"]["metrics_path"]).is_file()
assert all(Path(path).is_file() for path in qc["bwa_mem2_index"]["files"])
assert not Path(str(reference) + ".0123").exists()
PY
echo "[PASS] Stage 04 default BWA-MEM2 entrypoint and shared index cache"
