from __future__ import annotations

import json
import os
import shutil
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pysam
import zarr
from numcodecs import Blosc

from .cigar_encoding import pack_cigar, to_uint32
from .read_encoding import count_non_acgt, decode_read, pack_reads
from .timing import StepTimer


ENCODINGS = ("2bit_acgt", "4bit_acgtnp", "2bit_acgt_n_sparse")
READ_ALLOWED_SYMBOLS = frozenset({"A", "C", "G", "T", "N"})


@dataclass
class PrototypeConfig:
    sample_id: str
    alignment: Path
    contig: str
    output_dir: Path
    min_mapq: int = 30
    max_reads: int | None = None
    exclude_duplicate: bool = True
    exclude_supplementary: bool = False


@dataclass
class SegmentTables:
    read_id: list[int]
    chrom_id: list[int]
    ref_start: list[int]
    ref_end: list[int]
    query_start: list[int]
    query_end: list[int]
    strand: list[int]
    cigar_offset: list[int]
    cigar_len: list[int]
    cigar_ops: list[int]


def run_prototype(config: PrototypeConfig) -> dict[str, Any]:
    timer = StepTimer()
    config.output_dir.mkdir(parents=True, exist_ok=True)

    with timer.step("scan_alignment"):
        scan = scan_alignment(config)

    encoding_summaries = []
    for encoding in ENCODINGS:
        with timer.step(f"write_{encoding}"):
            encoding_summaries.append(write_encoding_store(
                output_dir=config.output_dir,
                sample_id=config.sample_id,
                contig=config.contig,
                encoding=encoding,
                scan=scan,
            ))

    with timer.step("write_summary"):
        summary = {
            "sample_id": config.sample_id,
            "contig": config.contig,
            "alignment": str(config.alignment),
            "filter": {
                "min_mapq": config.min_mapq,
                "exclude_duplicate": config.exclude_duplicate,
                "exclude_supplementary": config.exclude_supplementary,
            },
            "scan": scan["summary"],
            "encodings": encoding_summaries,
            "timing": timer.records,
        }
        write_json(config.output_dir / "summary.json", summary)
        write_tsv(config.output_dir / "summary.tsv", summary)
        write_json(config.output_dir / "timing.json", {"timing": timer.records})

    return summary


def scan_alignment(config: PrototypeConfig) -> dict[str, Any]:
    read_sequences: list[str] = []
    read_hp: list[int] = []
    read_ps: list[int] = []
    read_key_to_id: dict[tuple[str, int], int] = {}
    sequence_conflicts = 0
    hp_conflicts = 0
    ps_conflicts = 0

    tables = SegmentTables([], [], [], [], [], [], [], [], [], [])
    counts = Counter()
    hp_counts = Counter()
    ps_counts = Counter()
    total_bases = 0
    non_acgt_bases = 0
    n_bases = 0
    reads_with_n = 0
    reads_with_non_acgt = 0

    with pysam.AlignmentFile(str(config.alignment), "rb") as bam:
        bam_contig = resolve_contig(bam, config.contig)
        chrom_id = contig_to_id(config.contig)
        for read in bam.fetch(bam_contig):
            reason = rejection_reason(read, config)
            if reason:
                counts[f"excluded_{reason}"] += 1
                continue

            sequence = (read.query_sequence or "").upper()
            if not sequence:
                counts["excluded_missing_sequence"] += 1
                continue

            key = read_key(read)
            hp = read_haplotype(read)
            ps = read_phase_set(read)
            read_id = read_key_to_id.get(key)
            if read_id is None:
                read_id = len(read_sequences)
                read_key_to_id[key] = read_id
                read_sequences.append(sequence)
                read_hp.append(hp)
                read_ps.append(ps)
                total_bases += len(sequence)
                local_non_acgt, local_n = count_non_acgt(sequence)
                non_acgt_bases += local_non_acgt
                n_bases += local_n
                if local_n:
                    reads_with_n += 1
                if local_non_acgt:
                    reads_with_non_acgt += 1
            else:
                if read_sequences[read_id] != sequence:
                    sequence_conflicts += 1
                if read_hp[read_id] != hp:
                    read_hp[read_id] = 0
                    hp_conflicts += 1
                if read_ps[read_id] != ps:
                    read_ps[read_id] = 0
                    ps_conflicts += 1

            query_start = int(read.query_alignment_start or 0)
            query_end = int(read.query_alignment_end or len(sequence))
            ref_start = int(read.reference_start)
            ref_end = int(read.reference_end or read.reference_start)
            cigar_offset = len(tables.cigar_ops)
            cigar_ops = pack_cigar(read.cigartuples)
            tables.cigar_ops.extend(cigar_ops)

            tables.read_id.append(read_id)
            tables.chrom_id.append(chrom_id)
            tables.ref_start.append(ref_start)
            tables.ref_end.append(ref_end)
            tables.query_start.append(query_start)
            tables.query_end.append(query_end)
            tables.strand.append(-1 if read.is_reverse else 1)
            tables.cigar_offset.append(cigar_offset)
            tables.cigar_len.append(len(cigar_ops))

            hp_counts[str(hp)] += 1
            ps_counts["present" if ps else "missing"] += 1
            counts["kept_segments"] += 1
            if config.max_reads is not None and counts["kept_segments"] >= config.max_reads:
                break

    read_count = len(read_sequences)
    summary = {
        "read_count": read_count,
        "segment_count": len(tables.read_id),
        "total_read_bases": total_bases,
        "non_acgt_bases": non_acgt_bases,
        "n_bases": n_bases,
        "non_acgt_fraction": (non_acgt_bases / total_bases) if total_bases else 0.0,
        "n_fraction": (n_bases / total_bases) if total_bases else 0.0,
        "reads_with_n": reads_with_n,
        "reads_with_n_fraction": (reads_with_n / read_count) if read_count else 0.0,
        "reads_with_non_acgt": reads_with_non_acgt,
        "reads_with_non_acgt_fraction": (reads_with_non_acgt / read_count) if read_count else 0.0,
        "hp_segment_counts": dict(sorted(hp_counts.items())),
        "ps_segment_counts": dict(sorted(ps_counts.items())),
        "filter_counts": dict(sorted(counts.items())),
        "sequence_conflicts": sequence_conflicts,
        "hp_conflicts": hp_conflicts,
        "ps_conflicts": ps_conflicts,
    }
    print(
        "[SCAN] "
        f"reads={read_count} segments={len(tables.read_id)} bases={total_bases} "
        f"N={n_bases} ({summary['n_fraction']:.6%}) nonACGT={non_acgt_bases} "
        f"HP={summary['hp_segment_counts']} PS={summary['ps_segment_counts']}",
        flush=True,
    )
    return {
        "read_sequences": read_sequences,
        "read_hp": np.asarray(read_hp, dtype=np.uint8),
        "read_ps": np.asarray(read_ps, dtype=np.int64),
        "segments": tables,
        "summary": summary,
    }


def write_encoding_store(
    *,
    output_dir: Path,
    sample_id: str,
    contig: str,
    encoding: str,
    scan: dict[str, Any],
) -> dict[str, Any]:
    target = output_dir / encoding / "readset.zarr"
    if target.exists():
        shutil.rmtree(target)
    target.parent.mkdir(parents=True, exist_ok=True)

    allowed_symbols = READ_ALLOWED_SYMBOLS if encoding in {"4bit_acgtnp", "2bit_acgt_n_sparse"} else None
    packed = pack_reads(scan["read_sequences"], encoding, allowed_symbols=allowed_symbols)
    tables: SegmentTables = scan["segments"]
    compressor = Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE)
    root = zarr.open_group(str(target), mode="w")
    root.attrs.update({
        "schema": "mfm_integrated_readset_prototype",
        "schema_version": "prototype_2026_08_07",
        "sample_id": sample_id,
        "contig": contig,
        "encoding": encoding,
        **packed.attrs,
    })

    reads = root.create_group("reads")
    write_array(reads, "bases_packed", packed.bases, chunks=(max(1, min(len(packed.bases), 1 << 20)),), compressor=compressor)
    write_array(reads, "offsets", packed.offsets, chunks=(min(len(packed.offsets), 8192),), compressor=compressor)
    write_array(reads, "lengths", packed.lengths, chunks=(min(len(packed.lengths), 8192),), compressor=compressor)
    write_array(reads, "hp", scan["read_hp"], chunks=(min(len(scan["read_hp"]), 8192),), compressor=compressor)
    write_array(reads, "ps", scan["read_ps"], chunks=(min(len(scan["read_ps"]), 8192),), compressor=compressor)
    if packed.sparse_offsets is not None:
        write_array(reads, "sparse_offsets", packed.sparse_offsets, chunks=(min(len(packed.sparse_offsets), 8192),), compressor=compressor)
    if packed.sparse_positions is not None:
        write_array(reads, "sparse_positions", packed.sparse_positions, chunks=(max(1, min(len(packed.sparse_positions), 1 << 20)),), compressor=compressor)
    if packed.sparse_codes is not None:
        write_array(reads, "sparse_codes", packed.sparse_codes, chunks=(max(1, min(len(packed.sparse_codes), 1 << 20)),), compressor=compressor)

    segments = root.create_group("segments")
    write_array(segments, "read_id", np.asarray(tables.read_id, dtype=np.int32), chunks=(min(len(tables.read_id), 8192),), compressor=compressor)
    write_array(segments, "chrom_id", np.asarray(tables.chrom_id, dtype=np.int16), chunks=(min(len(tables.chrom_id), 8192),), compressor=compressor)
    write_array(segments, "ref_start", np.asarray(tables.ref_start, dtype=np.int32), chunks=(min(len(tables.ref_start), 8192),), compressor=compressor)
    write_array(segments, "ref_end", np.asarray(tables.ref_end, dtype=np.int32), chunks=(min(len(tables.ref_end), 8192),), compressor=compressor)
    write_array(segments, "query_start", np.asarray(tables.query_start, dtype=np.int32), chunks=(min(len(tables.query_start), 8192),), compressor=compressor)
    write_array(segments, "query_end", np.asarray(tables.query_end, dtype=np.int32), chunks=(min(len(tables.query_end), 8192),), compressor=compressor)
    write_array(segments, "strand", np.asarray(tables.strand, dtype=np.int8), chunks=(min(len(tables.strand), 8192),), compressor=compressor)
    write_array(segments, "cigar_offset", np.asarray(tables.cigar_offset, dtype=np.int64), chunks=(min(len(tables.cigar_offset), 8192),), compressor=compressor)
    write_array(segments, "cigar_len", np.asarray(tables.cigar_len, dtype=np.int32), chunks=(min(len(tables.cigar_len), 8192),), compressor=compressor)

    cigar = root.create_group("cigar")
    write_array(cigar, "ops_packed", to_uint32(tables.cigar_ops), chunks=(max(1, min(len(tables.cigar_ops), 1 << 20)),), compressor=compressor)

    decode_summary = validate_decoding(packed, scan["read_sequences"], encoding)
    size_bytes = directory_size(target)
    summary = {
        "encoding": encoding,
        "path": str(target),
        "size_bytes": size_bytes,
        "size_human": human_bytes(size_bytes),
        "read_count": len(scan["read_sequences"]),
        "segment_count": len(tables.read_id),
        "packed_attrs": packed.attrs,
        "decode": decode_summary,
    }
    write_json(target.parent / "encoding_summary.json", summary)
    print(
        f"[ENCODING] {encoding} size={summary['size_human']} "
        f"lossless={packed.attrs.get('lossless')} decode={decode_summary['status']}",
        flush=True,
    )
    return summary


def write_array(group: zarr.Group, name: str, data: np.ndarray, *, chunks: tuple[int, ...], compressor: Blosc) -> None:
    if data.size == 0:
        chunks = (1,)
    group.create_dataset(name, data=data, chunks=chunks, compressor=compressor, overwrite=True)


def validate_decoding(packed, sequences: list[str], encoding: str) -> dict[str, Any]:
    mismatch = 0
    checked = 0
    skipped_non_acgt = 0
    for i, sequence in enumerate(sequences[: min(len(sequences), 200)]):
        if encoding == "2bit_acgt" and count_non_acgt(sequence)[0] > 0:
            skipped_non_acgt += 1
            continue
        decoded = decode_read(packed, i, encoding)
        checked += 1
        if decoded != sequence:
            mismatch += 1
    return {
        "status": "pass" if mismatch == 0 else "fail",
        "checked_reads": checked,
        "mismatch_reads": mismatch,
        "skipped_non_acgt_reads": skipped_non_acgt,
        "note": "2bit_acgt is lossless only for A/C/G/T reads" if encoding == "2bit_acgt" else "",
    }


def rejection_reason(read: pysam.AlignedSegment, config: PrototypeConfig) -> str | None:
    if read.is_unmapped:
        return "unmapped"
    if read.is_secondary:
        return "secondary"
    if config.exclude_supplementary and read.is_supplementary:
        return "supplementary"
    if read.is_qcfail:
        return "qcfail"
    if config.exclude_duplicate and read.is_duplicate:
        return "duplicate"
    if int(read.mapping_quality) < config.min_mapq:
        return "low_mapq"
    if not read.query_sequence:
        return "missing_sequence"
    return None


def resolve_contig(bam: pysam.AlignmentFile, contig: str) -> str:
    references = set(bam.references)
    if contig in references:
        return contig
    alt = contig.removeprefix("chr")
    if alt in references:
        return alt
    alt = f"chr{contig}"
    if alt in references:
        return alt
    raise ValueError(f"contig not found in alignment: {contig}")


def contig_to_id(contig: str) -> int:
    value = contig.removeprefix("chr")
    if value == "X":
        return 23
    if value == "Y":
        return 24
    if value in {"M", "MT"}:
        return 25
    return int(value)


def mate_code(read: pysam.AlignedSegment) -> int:
    return 1 if read.is_read1 else (2 if read.is_read2 else 0)


def read_key(read: pysam.AlignedSegment) -> tuple[str, int]:
    return (str(read.query_name), mate_code(read))


def read_haplotype(read: pysam.AlignedSegment) -> int:
    try:
        value = int(read.get_tag("HP"))
    except (KeyError, TypeError, ValueError):
        return 0
    return value if value in (1, 2) else 0


def read_phase_set(read: pysam.AlignedSegment) -> int:
    try:
        value = int(read.get_tag("PS"))
    except (KeyError, TypeError, ValueError):
        return 0
    return value if value > 0 else 0


def directory_size(path: Path) -> int:
    total = 0
    for root, _dirs, files in os.walk(path):
        for file_name in files:
            total += (Path(root) / file_name).stat().st_size
    return total


def human_bytes(size: int) -> str:
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return f"{value:.2f} {unit}"
        value /= 1024.0
    return f"{size} B"


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")


def write_tsv(path: Path, summary: dict[str, Any]) -> None:
    rows = []
    scan = summary["scan"]
    for enc in summary["encodings"]:
        rows.append({
            "sample_id": summary["sample_id"],
            "contig": summary["contig"],
            "encoding": enc["encoding"],
            "size_bytes": enc["size_bytes"],
            "size_human": enc["size_human"],
            "read_count": scan["read_count"],
            "segment_count": scan["segment_count"],
            "total_read_bases": scan["total_read_bases"],
            "n_bases": scan["n_bases"],
            "n_fraction": scan["n_fraction"],
            "non_acgt_bases": scan["non_acgt_bases"],
            "non_acgt_fraction": scan["non_acgt_fraction"],
            "lossless": enc["packed_attrs"].get("lossless"),
            "decode_status": enc["decode"]["status"],
        })
    columns = list(rows[0]) if rows else []
    with path.open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[col]) for col in columns) + "\n")
