from __future__ import annotations

import csv
import gzip
import json
import os
import shlex
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


DEFAULT_MINIMAP2 = "minimap2"


@dataclass
class AssemblyAnchorConfig:
    cohort: str
    output_root: Path
    reference_fasta: Path
    sample: str
    assembly_dir: Path | None = None
    label: str = ""
    minimap2_bin: str = DEFAULT_MINIMAP2
    preset: str = "asm20"
    threads: int = 48
    min_mapq: int = 30
    min_contig_len: int = 100_000
    overwrite: bool = False
    dry_run: bool = False
    summary_tsv: Path | None = None
    metadata_json: Path | None = None


@dataclass
class AnchorRow:
    sample: str
    label: str
    fasta: str
    paf: str
    total_contigs: int
    total_bases: int
    anchored_contigs: int
    anchored_query_bases: int
    anchored_fraction: float
    best_target_count: int
    elapsed_sec: float
    status: str
    message: str = ""


def run_assembly_anchor(config: AssemblyAnchorConfig) -> dict[str, Any]:
    started = time.time()
    validate_config(config)
    assembly_dir = resolve_assembly_dir(config)
    fasta_dir = assembly_dir / "fastas"
    anchor_dir = assembly_dir / "anchor"
    anchor_dir.mkdir(parents=True, exist_ok=True)
    fastas = sorted(fasta_dir.glob("*.fa")) + sorted(fasta_dir.glob("*.fa.gz"))
    if not fastas:
        raise FileNotFoundError(f"no assembly FASTA found under {fasta_dir}")

    rows: list[AnchorRow] = []
    for fasta in fastas:
        rows.append(anchor_one(config, assembly_dir, anchor_dir, fasta))

    summary_tsv = config.summary_tsv or (anchor_dir / "assembly_anchor.summary.tsv")
    metadata_json = config.metadata_json or (anchor_dir / "assembly_anchor.metadata.json")
    write_summary(summary_tsv, rows)
    metadata = {
        "schema": "mfm_assembly_anchor_v1",
        "cohort": config.cohort,
        "sample": config.sample,
        "label": config.label,
        "assembly_dir": str(assembly_dir),
        "reference_fasta": str(config.reference_fasta),
        "preset": config.preset,
        "min_mapq": config.min_mapq,
        "min_contig_len": config.min_contig_len,
        "threads": config.threads,
        "minimap2_bin": config.minimap2_bin,
        "started_at_unix": started,
        "completed_at_unix": time.time(),
        "elapsed_sec": round(time.time() - started, 3),
        "summary_tsv": str(summary_tsv),
        "ok_count": sum(1 for row in rows if row.status in {"OK", "SKIP"}),
        "fail_count": sum(1 for row in rows if row.status == "FAIL"),
    }
    write_json(metadata_json, metadata)
    if metadata["fail_count"]:
        raise RuntimeError(f"assembly anchoring failed; see {summary_tsv}")
    return metadata


def anchor_one(config: AssemblyAnchorConfig, assembly_dir: Path, anchor_dir: Path, fasta: Path) -> AnchorRow:
    t0 = time.time()
    paf = anchor_dir / f"{fasta.name}.to_reference.paf"
    try:
        if config.dry_run:
            paf.write_text("", encoding="utf-8")
            return make_row(config, fasta, paf, t0, "OK", "dry-run")
        if not paf.exists() or config.overwrite:
            cmd = [config.minimap2_bin, "-x", config.preset, "-t", str(max(1, config.threads)), str(config.reference_fasta), str(fasta)]
            run_to_file(cmd, paf, anchor_dir / f"{fasta.name}.minimap2.log")
        return make_row(config, fasta, paf, t0, "OK", "ok")
    except Exception as exc:  # noqa: BLE001
        return AnchorRow(
            sample=config.sample,
            label=config.label,
            fasta=str(fasta),
            paf=str(paf),
            total_contigs=0,
            total_bases=0,
            anchored_contigs=0,
            anchored_query_bases=0,
            anchored_fraction=0.0,
            best_target_count=0,
            elapsed_sec=round(time.time() - t0, 3),
            status="FAIL",
            message=f"{type(exc).__name__}: {exc}",
        )


def make_row(config: AssemblyAnchorConfig, fasta: Path, paf: Path, t0: float, status: str, message: str) -> AnchorRow:
    lengths = fasta_lengths(fasta)
    total_bases = sum(lengths.values())
    best = parse_best_paf(paf, config.min_mapq) if paf.exists() else {}
    anchored_bases = sum(end - start for start, end, _target, _mapq in best.values())
    targets = {target for _start, _end, target, _mapq in best.values()}
    return AnchorRow(
        sample=config.sample,
        label=config.label,
        fasta=str(fasta),
        paf=str(paf),
        total_contigs=len(lengths),
        total_bases=total_bases,
        anchored_contigs=len(best),
        anchored_query_bases=anchored_bases,
        anchored_fraction=round(anchored_bases / total_bases, 6) if total_bases else 0.0,
        best_target_count=len(targets),
        elapsed_sec=round(time.time() - t0, 3),
        status=status,
        message=message,
    )


def parse_best_paf(paf: Path, min_mapq: int) -> dict[str, tuple[int, int, str, int]]:
    best: dict[str, tuple[int, int, str, int]] = {}
    best_score: dict[str, tuple[int, int]] = {}
    with paf.open("rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname = fields[0]
            qstart = int(fields[2])
            qend = int(fields[3])
            target = fields[5]
            mapq = int(fields[11])
            if mapq < min_mapq:
                continue
            aln = max(0, qend - qstart)
            score = (aln, mapq)
            if score > best_score.get(qname, (-1, -1)):
                best_score[qname] = score
                best[qname] = (qstart, qend, target, mapq)
    return best


def fasta_lengths(fasta: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name: str | None = None
    current = 0
    opener = gzip.open if fasta.suffix == ".gz" else open
    with opener(fasta, "rt") as fh:  # type: ignore[arg-type]
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = current
                name = line[1:].strip().split()[0]
                current = 0
            else:
                current += len(line.strip())
    if name is not None:
        lengths[name] = current
    return lengths


def validate_config(config: AssemblyAnchorConfig) -> None:
    if not config.reference_fasta.exists():
        raise FileNotFoundError(f"reference FASTA not found: {config.reference_fasta}")
    path = Path(config.minimap2_bin)
    if not ((path.exists() and os.access(path, os.X_OK)) or shutil_which(config.minimap2_bin)):
        raise FileNotFoundError(f"minimap2 not found: {config.minimap2_bin}")


def resolve_assembly_dir(config: AssemblyAnchorConfig) -> Path:
    if config.assembly_dir is not None:
        path = config.assembly_dir
    else:
        if not config.label:
            raise ValueError("provide --assembly-dir or --label")
        path = config.output_root / config.sample / "denovo_assembly" / config.label
    if not path.exists():
        raise FileNotFoundError(f"assembly dir not found: {path}")
    return path.resolve()


def run_to_file(cmd: list[str], stdout_path: Path, log_path: Path) -> None:
    env = os.environ.copy()
    exe_parent = str(Path(cmd[0]).resolve().parent) if "/" in cmd[0] else ""
    if exe_parent:
        env["PATH"] = exe_parent + os.pathsep + env.get("PATH", "")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("wt", encoding="utf-8") as out, log_path.open("wt", encoding="utf-8") as log:
        log.write("$ " + " ".join(shlex.quote(item) for item in cmd) + "\n")
        log.flush()
        proc = subprocess.run(cmd, stdout=out, stderr=log, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed with exit code {proc.returncode}; see {log_path}")


def write_summary(path: Path, rows: list[AnchorRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(asdict(rows[0]).keys()) if rows else [field for field in AnchorRow.__dataclass_fields__]
    with path.open("wt", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def shutil_which(value: str) -> str | None:
    from shutil import which

    return which(value)
