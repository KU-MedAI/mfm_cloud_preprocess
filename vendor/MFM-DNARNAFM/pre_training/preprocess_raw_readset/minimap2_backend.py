from __future__ import annotations

import subprocess
from pathlib import Path


def run_minimap2(reference_fasta: Path, individual_fasta: Path, paf: Path,
                 executable: str = "minimap2") -> None:
    """Create an assembly-to-reference PAF with CIGAR for non-VCF-derived FASTA."""
    paf.parent.mkdir(parents=True, exist_ok=True)
    with paf.open("w") as out:
        result = subprocess.run(
            [
                executable,
                "-c",
                "-x",
                "asm5",
                "--secondary=no",
                str(reference_fasta),
                str(individual_fasta),
            ],
            stdout=out,
            stderr=subprocess.PIPE,
            text=True,
        )
    if result.returncode or not paf.stat().st_size:
        raise RuntimeError(f"minimap2 failed: {result.stderr}")
