"""Record TSV helpers for shared-HAP preprocessing.

The preprocessing CLI consumes one TSV row per sample/contig. Keeping TSV
validation here makes wrapper-generated and manually-authored manifests follow
the same rules.
"""

from __future__ import annotations

import csv
from pathlib import Path

RECORD_COLUMNS = (
    "sample",
    "contig",
    "hap1_fasta",
    "hap1_chain",
    "hap2_fasta",
    "hap2_chain",
    "bam",
)


def read_record_tsv(path: Path) -> list[dict]:
    """Read and validate a shared-HAP record TSV."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="	")
        missing = set(RECORD_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"record TSV missing columns: {sorted(missing)}")
        records = list(reader)
    if not records:
        raise ValueError(f"record TSV has no rows: {path}")
    return records


def hap_paths(job: dict, label: str) -> tuple[Path, Path | None]:
    """Return FASTA/chain paths for HAP1 or HAP2 from a job dict.

    SNP-only cohorts preserve reference coordinates and intentionally do not
    produce chain files. An empty chain column therefore means identity
    coordinate mapping, not a missing required file.
    """
    key = label.lower()
    chain_value = str(job.get(f"{key}_chain", "") or "")
    return Path(job[f"{key}_fasta"]), (Path(chain_value) if chain_value else None)
