"""Shared-read, haplotype-specific coordinate zarr writer.

This module implements the schema chosen for read-base-aware model work:

* one chromosome zarr per sample/contig;
* one shared read table per chromosome, so read bases are not duplicated across
  HAP1 and HAP2;
* one HAP1 group and one HAP2 group, each with its own sequence windows and
  precomputed window-local read segment coordinates.

The dataloader should only need to read a haplotype sequence, read indices, and
query_start/query_end coordinates, then slice the shared read sequence. Expensive
BAM CIGAR traversal and GRCh38->haplotype coordinate conversion are completed in
preprocessing.
"""

from __future__ import annotations

import os

import bisect
import csv
import json
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pysam
import zarr
from numcodecs import Blosc, VLenUTF8

from .config import ExportConfig
from .coordinates import (
    ChainAlignment,
    ProjectedWindowSegment,
    chain_alignment_to_map,
    identity_coordinate_map,
    interval_mapping_coverage,
    load_chain_alignment,
    project_read_to_window_segment,
    project_read_to_window_segment_blockwise,
    reference_interval,
)
from .exporter import _fasta_sequence
from .filters import rejection_reason
from .shared_hap.manifest import utc_now, write_json
from .shared_hap.schema import (
    COMPATIBLE_SHARED_SCHEMAS,
    HAP_LABELS,
    READ_CHUNK_SIZE,
    READ_HAPLOTYPE_STORAGE,
    READ_PHASE_SET_STORAGE,
    READ_SEQUENCE_STORAGE,
    READ_CIGAR_STORAGE,
    READ_REFERENCE_STORAGE,
    SCHEMA_VERSION,
    SEGMENT_CHUNK_SIZE,
    SHARED_SCHEMA,
    VALIDATION_LEVELS,
)
from .windowing import deterministic_offset, make_windows, trim_terminal_n
from .window_filter_stats import DROP_REASON_PRIORITY, count_index_exclusions, index_filter_policy




def _timing_add(timing: dict[str, Any] | None, name: str, elapsed: float) -> None:
    if timing is None:
        return
    steps = timing.setdefault("steps", {})
    steps[name] = float(steps.get(name, 0.0)) + float(elapsed)


def _timing_set_count(timing: dict[str, Any] | None, name: str, value: int) -> None:
    if timing is None:
        return
    counts = timing.setdefault("counts", {})
    counts[name] = int(value)


HAP_PROVENANCE_ATTR_KEYS = (
    "coordinate_system",
    "coordinate_map_source",
    "source_chain_path",
    "source_individual_fasta",
    "source_bam_path",
    "window_index_path",
    "window_index_schema",
    "window_index_kept_count",
    "window_index_total_count",
    "window_index_dropped_n",
    "window_index_dropped_blacklist",
    "window_index_drop_reason_priority",
    "random_offset",
    "random_offset_unit",
    "sequence_valid_window_count_total",
    "coordinate_eligible_window_count_total",
    "coordinate_excluded_window_count",
    "coordinate_excluded_bases",
)


def _mate_code(read: pysam.AlignedSegment) -> int:
    return 1 if read.is_read1 else (2 if read.is_read2 else 0)


def _read_key(read: pysam.AlignedSegment) -> tuple[str, int]:
    return (str(read.query_name), _mate_code(read))


def _read_haplotype(read: pysam.AlignedSegment) -> int:
    """Return HP tag as 0/1/2, using 0 for untagged or non-standard values."""
    try:
        value = int(read.get_tag("HP"))
    except (KeyError, TypeError, ValueError):
        return 0
    return value if value in (1, 2) else 0


def _read_phase_set(read: pysam.AlignedSegment) -> int:
    """Return the BAM PS tag as a phase-set id, using 0 for missing/invalid.

    In HP/PS-tagged BAMs, HP indicates the haplotype assignment and PS indicates
    the phase set/block that supports the assignment. PS is usually represented
    as a genomic coordinate-like integer. We keep the raw integer because the
    zarr file is already chromosome-scoped; 0 is reserved for unknown.
    """
    try:
        value = int(read.get_tag("PS"))
    except (KeyError, TypeError, ValueError):
        return 0
    return value if value > 0 else 0


def _alignment_file(bam_path: Path, reference_fasta: Path | None = None) -> pysam.AlignmentFile:
    kwargs = {}
    if reference_fasta is not None:
        kwargs["reference_filename"] = str(reference_fasta)
    threads=int(os.environ.get('MFM_COMPACT_ALIGNMENT_THREADS','1'))
    if threads<1:raise ValueError('MFM_COMPACT_ALIGNMENT_THREADS는 양수여야 함')
    return pysam.AlignmentFile(str(bam_path), "rb", threads=threads, **kwargs)


def _clipping_removed_query_sequence(read: pysam.AlignedSegment) -> tuple[str, int, int]:
    """Return mapped query sequence and its original query-coordinate bounds.

    pysam's query_alignment_start/end exclude soft clips. Hard clips are absent
    from query_sequence. The returned sequence preserves M/=/X and insertion
    bases while removing clipping bases.
    """
    sequence = (read.query_sequence or "").upper()
    start = int(read.query_alignment_start or 0)
    end = int(read.query_alignment_end or len(sequence))
    return sequence[start:end], start, end


def _pack_cigartuples(cigartuples: list[tuple[int, int]] | tuple[tuple[int, int], ...] | None) -> list[int]:
    """Pack CIGAR as uint32 values matching the stored clipping-removed query.

    HTSlib encodes one CIGAR op as (length << 4) | op_code. The stored read
    sequence in this writer has soft/hard clipping removed, so clipping ops are
    omitted here as well. This keeps query_start/query_end and CIGAR operations
    in the same query coordinate system used by the zarr read_sequence.
    """
    if not cigartuples:
        return []
    packed: list[int] = []
    for op_code, length in cigartuples:
        op_code = int(op_code)
        length = int(length)
        if op_code in (4, 5):  # soft/hard clipping, not present in stored query sequence
            continue
        if length < 0:
            raise ValueError(f"negative CIGAR length: {length}")
        if length >= (1 << 28):
            raise ValueError(f"CIGAR length too large for uint32 packing: {length}")
        packed.append((length << 4) | op_code)
    return packed


def summarize_hap_windows(
    *,
    sample: str,
    label: str,
    contig: str,
    individual_fasta: Path,
    chain_path: Path | None,
    config: ExportConfig,
) -> dict:
    """Summarize valid/N-filtered and coordinate-eligible windows for one hap."""
    config.validate()
    query_name, original = _fasta_sequence(individual_fasta)
    trimmed = trim_terminal_n(original)
    offset_unit = config.effective_offset_unit()
    offset = deterministic_offset(config.seed, sample, f"{label}:{contig}", offset_unit)
    windows = make_windows(
        trimmed.sequence,
        source_start=trimmed.source_start,
        offset=offset,
        window_size=config.window_size,
        max_n_ratio=config.max_n_ratio,
    )
    chain_alignment = _load_optional_chain_alignment(
        chain_path,
        query_name=query_name,
        target_name=contig,
        raw_length=len(original),
    )
    mapped_intervals = _mapped_intervals_from_alignment(chain_alignment, len(original))
    rows = [
        (*interval_mapping_coverage(mapped_intervals, window.raw_start, window.raw_end), window)
        for window in windows
    ]
    eligible = [row for row in rows if row[1] >= config.min_coordinate_coverage]
    coverage = [float(row[1]) for row in rows]
    return {
        "sample": sample,
        "label": label,
        "contig": contig,
        "query_name": query_name,
        "terminal_n_left": trimmed.source_start,
        "terminal_n_right": len(original) - trimmed.source_end,
        "random_offset": offset,
        "random_offset_unit": offset_unit,
        "sequence_valid_window_count_total": len(windows),
        "coordinate_eligible_window_count_total": len(eligible),
        "coordinate_excluded_window_count": len(windows) - len(eligible),
        "coordinate_excluded_bases": sum(
            row[2].actual_length for row in rows
            if row[1] < config.min_coordinate_coverage
        ),
        "coordinate_coverage_min": min(coverage) if coverage else None,
        "coordinate_coverage_max": max(coverage) if coverage else None,
        "coordinate_coverage_threshold": config.min_coordinate_coverage,
        "coordinate_map_source": _coordinate_map_source(chain_alignment),
        "source_chain_path": str(chain_path) if chain_path else "",
    }


def build_hap_window_artifacts(
    *,
    sample: str,
    label: str,
    contig: str,
    individual_fasta: Path,
    bam_path: Path,
    chain_path: Path | None,
    reference_fasta: Path | None = None,
    max_windows: int | None = None,
    window_start_index: int = 0,
    config: ExportConfig = ExportConfig(),
) -> dict:
    """Build one HAP chunk using shared-read references instead of read strings.

    Returned records contain enough information for the final zarr writer to
    register each clipping-removed read sequence once, then store only
    window-local segment coordinates under the HAP group.
    """
    config.validate()
    if window_start_index < 0:
        raise ValueError("window_start_index must be non-negative")
    query_name, original = _fasta_sequence(individual_fasta)
    trimmed = trim_terminal_n(original)
    offset_unit = config.effective_offset_unit()
    offset = deterministic_offset(config.seed, sample, f"{label}:{contig}", offset_unit)
    all_windows = make_windows(
        trimmed.sequence,
        source_start=trimmed.source_start,
        offset=offset,
        window_size=config.window_size,
        max_n_ratio=config.max_n_ratio,
    )
    chain_alignment = _load_optional_chain_alignment(
        chain_path,
        query_name=query_name,
        target_name=contig,
        raw_length=len(original),
    )
    mapped_intervals = _mapped_intervals_from_alignment(chain_alignment, len(original))
    coverage_rows = [
        (*interval_mapping_coverage(mapped_intervals, window.raw_start, window.raw_end), window)
        for window in all_windows
    ]
    eligible_rows = [row for row in coverage_rows if row[1] >= config.min_coordinate_coverage]
    stop = None if max_windows is None else window_start_index + max_windows
    selected_rows = eligible_rows[window_start_index:stop]
    windows = [row[2] for row in selected_rows]
    if not windows:
        raise ValueError(
            f"no coordinate-eligible window selected: sample={sample} label={label} "
            f"contig={contig} start={window_start_index} eligible={len(eligible_rows)}"
        )

    flank = 1_000_000
    raw_range = (
        max(0, windows[0].raw_start - flank),
        min(len(original), windows[-1].raw_end + flank),
    )
    if _coordinate_map_source(chain_alignment) == "identity":
        coord = identity_coordinate_map(len(original), raw_range=raw_range)
    else:
        coord = chain_alignment_to_map(chain_alignment, len(original), raw_range=raw_range)

    records: list[list[dict[str, Any]]] = [[] for _ in windows]
    seen_by_window: list[set[tuple[Any, ...]]] = [set() for _ in windows]
    counts = Counter()

    with _alignment_file(bam_path, reference_fasta) as bam:
        bam_contig = contig if contig in bam.references else contig.removeprefix("chr")
        for wi, window in enumerate(windows):
            r0, r1 = reference_interval(coord, window.raw_start, window.raw_end)
            for read in bam.fetch(bam_contig, r0, r1):
                reason = rejection_reason(read, config)
                if reason:
                    counts[f"excluded_{reason}"] += 1
                    continue
                projected = project_read_to_window_segment(
                    read,
                    coord,
                    window.raw_start,
                    window.raw_start + window.actual_length,
                )
                if projected is None:
                    counts["excluded_unprojectable"] += 1
                    continue
                if projected.raw_start >= window.raw_end or projected.raw_end <= window.raw_start:
                    continue
                read_sequence, aln_query_start, aln_query_end = _clipping_removed_query_sequence(read)
                if not read_sequence:
                    counts["excluded_empty_mapped_query"] += 1
                    continue
                query_start = projected.query_start - aln_query_start
                query_end = projected.query_end - aln_query_start
                if query_start < 0 or query_end <= query_start or query_end > len(read_sequence):
                    counts["excluded_invalid_query_slice"] += 1
                    continue
                read_start = projected.window_start - window.raw_start
                read_end = projected.window_end - window.raw_start
                if read_start < 0 or read_end <= read_start or read_end > int(window.actual_length):
                    counts["excluded_invalid_hap_interval"] += 1
                    continue

                key = _read_key(read)
                shared_read_key = (
                    key,
                    int(read.reference_start),
                    int(read.reference_end or read.reference_start),
                    read.cigarstring or "",
                    bool(read.is_supplementary),
                    bool(read.is_reverse),
                )
                segment_key = (
                    shared_read_key,
                    int(query_start),
                    int(query_end),
                    int(read_start),
                    int(read_end),
                )
                if segment_key in seen_by_window[wi]:
                    continue
                seen_by_window[wi].add(segment_key)
                records[wi].append({
                    "read_key": key,
                    "shared_read_key": shared_read_key,
                    "read_sequence": read_sequence,
                    "read_haplotype": _read_haplotype(read),
                    "read_phase_set": _read_phase_set(read),
                    "read_index": -1,  # filled by SharedHapChrWriter
                    "query_start": int(query_start),
                    "query_end": int(query_end),
                    "read_start": int(read_start),
                    "read_end": int(read_end),
                    "read_strand": -1 if read.is_reverse else 1,
                    "reference_start": int(read.reference_start),
                    "reference_end": int(read.reference_end or read.reference_start),
                    "cigar": read.cigarstring or "",
                    "cigar_ops_packed": _pack_cigartuples(read.cigartuples),
                    "segment_type": 1 if read.is_supplementary else 0,
                })
                counts["kept_segments"] += 1

    attrs = {
        "sample": sample,
        "label": label,
        "contig": contig,
        "coordinate_system": "individual_raw_0based_halfopen",
        "terminal_n_left": trimmed.source_start,
        "terminal_n_right": len(original) - trimmed.source_end,
        "random_offset": offset,
        "random_offset_unit": offset_unit,
        "valid_window_start_index": window_start_index,
        "sequence_valid_window_count_total": len(all_windows),
        "coordinate_eligible_window_count_total": len(eligible_rows),
        "coordinate_excluded_window_count": len(all_windows) - len(eligible_rows),
        "coordinate_excluded_bases": sum(
            row[2].actual_length for row in coverage_rows
            if row[1] < config.min_coordinate_coverage
        ),
        "selected_coordinate_mapped_bases": [row[0] for row in selected_rows],
        "selected_coordinate_coverage": [row[1] for row in selected_rows],
        "coordinate_map_source": _coordinate_map_source(chain_alignment),
        "source_chain_path": str(chain_path) if chain_path else "",
        "source_individual_fasta": str(individual_fasta),
        "source_bam_path": str(bam_path),
        "config": config.to_dict(),
    }
    return {
        "label": label,
        "windows": windows,
        "records": records,
        "counts": counts,
        "attrs": attrs,
        "window_count": len(windows),
        "segment_count": sum(len(x) for x in records),
    }



def prepare_hap_context(
    *,
    sample: str,
    label: str,
    contig: str,
    individual_fasta: Path,
    chain_path: Path | None,
    config: ExportConfig,
    window_index_path: Path | None = None,
) -> dict[str, Any]:
    """Precompute per-HAP FASTA windows and coordinate coverage once per contig.

    The earlier single-HAP builder recomputed FASTA loading, terminal-N trimming,
    deterministic windows, chain mapped intervals, and coordinate coverage for
    every computation chunk. This context object keeps those read-only values in
    memory inside the contig worker and lets paired chunk processing reuse them.
    """
    config.validate()
    query_name, original = _fasta_sequence(individual_fasta)
    trimmed = trim_terminal_n(original)
    offset_unit = config.effective_offset_unit()
    offset = deterministic_offset(config.seed, sample, f"{label}:{contig}", offset_unit)
    all_windows = make_windows(
        trimmed.sequence,
        source_start=trimmed.source_start,
        offset=offset,
        window_size=config.window_size,
        max_n_ratio=config.max_n_ratio,
    )
    chain_alignment = _load_optional_chain_alignment(
        chain_path,
        query_name=query_name,
        target_name=contig,
        raw_length=len(original),
    )
    mapped_intervals = _mapped_intervals_from_alignment(chain_alignment, len(original))
    coverage_rows = [
        (*interval_mapping_coverage(mapped_intervals, window.raw_start, window.raw_end), window)
        for window in all_windows
    ]
    coordinate_eligible_rows = [row for row in coverage_rows if row[1] >= config.min_coordinate_coverage]
    eligible_rows, window_index_summary = _apply_external_window_index(
        window_index_path=window_index_path,
        label=label,
        contig=contig,
        coverage_rows=coordinate_eligible_rows,
    )
    coverage = [float(row[1]) for row in coverage_rows]
    summary = {
        "sample": sample,
        "label": label,
        "contig": contig,
        "query_name": query_name,
        "terminal_n_left": trimmed.source_start,
        "terminal_n_right": len(original) - trimmed.source_end,
        "random_offset": offset,
        "random_offset_unit": offset_unit,
        "sequence_valid_window_count_total": len(all_windows),
        "coordinate_eligible_window_count_total": len(eligible_rows),
        "coordinate_candidate_window_count_total": len(coordinate_eligible_rows),
        "coordinate_excluded_window_count": len(all_windows) - len(coordinate_eligible_rows),
        "coordinate_excluded_bases": sum(
            row[2].actual_length for row in coverage_rows
            if row[1] < config.min_coordinate_coverage
        ),
        **window_index_summary,
        "coordinate_coverage_min": min(coverage) if coverage else None,
        "coordinate_coverage_max": max(coverage) if coverage else None,
        "coordinate_coverage_threshold": config.min_coordinate_coverage,
        "coordinate_map_source": _coordinate_map_source(chain_alignment),
        "source_chain_path": str(chain_path) if chain_path else "",
    }
    return {
        "sample": sample,
        "label": label,
        "contig": contig,
        "individual_fasta": Path(individual_fasta),
        "chain_path": Path(chain_path) if chain_path else None,
        "chain_alignment": chain_alignment,
        "query_name": query_name,
        "original": original,
        "trimmed": trimmed,
        "all_windows": all_windows,
        "coverage_rows": coverage_rows,
        "eligible_rows": eligible_rows,
        "summary": summary,
    }



def _apply_external_window_index(
    *,
    window_index_path: Path | None,
    label: str,
    contig: str,
    coverage_rows: list[tuple[int, float, Any]],
) -> tuple[list[tuple[int, float, Any]], dict[str, Any]]:
    """Apply a precomputed 05 windows.tsv file to HAP context rows.

    The direct compact path still computes sequence windows from FASTA so the
    sequence payload remains source-of-truth.  When a 05 window index is passed,
    this function only replaces the eligible-window list with the `kept=1` rows
    from windows.tsv.  This is how process_code's blacklist/N window policy is
    propagated into the final model-input zarr without changing the zarr layout.
    """
    if window_index_path is None or not str(window_index_path):
        return coverage_rows, {
            "window_index_path": "",
            "window_index_schema": "",
            "window_index_total_count": 0,
            "window_index_kept_count": 0,
            "window_index_dropped_n": 0,
            "window_index_dropped_blacklist": 0,
            "window_index_drop_reason_priority": "",
        }
    path = Path(window_index_path)
    if path.is_dir():
        path = path / "windows.tsv"
    if not path.exists() or not path.stat().st_size:
        raise FileNotFoundError(f"window index TSV is missing: {path}")

    matching: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"contig", "haplotype", "raw_start", "actual_length", "kept"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"window index TSV missing columns: {sorted(missing)}: {path}")
        for row in reader:
            if row.get("contig") == contig and row.get("haplotype") == label:
                matching.append(dict(row))
    if not matching:
        raise ValueError(f"window index has no rows for {label}/{contig}: {path}")

    by_key: dict[tuple[int, int], tuple[int, float, Any]] = {}
    for row in coverage_rows:
        window = row[2]
        by_key[(int(window.raw_start), int(window.actual_length))] = row

    selected: list[tuple[int, float, Any]] = []
    missing_keys: list[tuple[int, int]] = []
    drop_counts = count_index_exclusions(matching, **index_filter_policy(path, matching))
    kept_count = 0
    for row in matching:
        kept = str(row.get("kept", "")).strip() in {"1", "true", "True", "TRUE"}
        if not kept:
            continue
        kept_count += 1
        key = (int(row["raw_start"]), int(row["actual_length"]))
        found = by_key.get(key)
        if found is None:
            missing_keys.append(key)
            continue
        selected.append(found)
    if missing_keys:
        preview = ", ".join(f"{start}:{length}" for start, length in missing_keys[:5])
        raise ValueError(
            f"window index kept rows do not match FASTA-derived windows for {label}/{contig}: {preview}. "
            "Check window_size/offset_unit/seed/drop_n consistency."
        )
    if kept_count <= 0:
        raise ValueError(f"window index kept no windows for {label}/{contig}: {path}")
    return selected, {
        "window_index_path": str(path),
        "window_index_schema": json.loads(path.with_name("windows.json").read_text()).get("schema", "mfm_integrated_window_index_v1") if path.with_name("windows.json").exists() else "mfm_integrated_window_index_v1",
        "window_index_total_count": len(matching),
        "window_index_kept_count": kept_count,
        "window_index_dropped_n": drop_counts["dropped_n"],
        "window_index_dropped_blacklist": drop_counts["dropped_blacklist"],
        "window_index_drop_reason_priority": DROP_REASON_PRIORITY,
    }

def _identity_chain_alignment(query_name: str, target_name: str, raw_length: int) -> ChainAlignment:
    """Return a synthetic full-length plus/plus alignment for SNP-only FASTA."""
    header = {
        "score": raw_length,
        "target_name": target_name,
        "target_size": raw_length,
        "target_strand": "+",
        "target_start": 0,
        "target_end": raw_length,
        "query_name": query_name,
        "query_size": raw_length,
        "query_strand": "+",
        "query_start": 0,
        "query_end": raw_length,
        "chain_id": "identity",
        "source": "identity",
    }
    return ChainAlignment(header=header, blocks=((0, 0, raw_length),))


def _load_optional_chain_alignment(
    chain_path: Path | None,
    *,
    query_name: str,
    target_name: str,
    raw_length: int,
) -> ChainAlignment:
    """Load a real chain, or use identity mapping when chain is intentionally absent.

    process_code uses SNP-only personal FASTA for GEUVADIS. SNP-only consensus
    preserves reference coordinates, so no chain file is generated or required.
    Empty chain paths are therefore treated as identity maps.
    """
    if chain_path is None or not str(chain_path):
        return _identity_chain_alignment(query_name, target_name, raw_length)
    if not chain_path.exists() or not chain_path.stat().st_size:
        raise FileNotFoundError(f"chain file is required unless SNP-only identity mapping is intended: {chain_path}")
    return load_chain_alignment(chain_path, query_name=query_name, target_name=target_name)


def _coordinate_map_source(alignment: ChainAlignment) -> str:
    return "identity" if alignment.header.get("source") == "identity" else "chain"


def _mapped_intervals_from_alignment(alignment: ChainAlignment, raw_length: int) -> list[tuple[int, int]]:
    if _coordinate_map_source(alignment) == "identity":
        return [(0, raw_length)]
    return _merge_intervals([(q, q + length) for q, _t, length in alignment.blocks])


def _coordinate_map_from_context(context: dict[str, Any], raw_range: tuple[int, int]):
    alignment = context["chain_alignment"]
    original = context["original"]
    if _coordinate_map_source(alignment) == "identity":
        return identity_coordinate_map(len(original), raw_range=raw_range)
    return chain_alignment_to_map(alignment, len(original), raw_range=raw_range)


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    intervals = sorted((int(s), int(e)) for s, e in intervals if int(e) > int(s))
    if not intervals:
        return []
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def _candidate_window_indices(
    *,
    interval_starts: list[int],
    interval_ends: list[int],
    read_ref_start: int,
    read_ref_end: int,
) -> range:
    """Return the sorted window-index range whose ref intervals overlap a read.

    Haplotype windows are generated in coordinate order, and the bcftools
    consensus chain is plus/plus, so their projected reference intervals remain
    monotonic. A long read normally overlaps only one or two 1M windows. Using
    binary search avoids re-projecting the same read against every window in a
    chunk, which was the dominant Python-side cost when `windows_per_chunk` was
    increased for better zarr write efficiency.
    """
    if read_ref_end <= read_ref_start:
        return range(0, 0)
    first = bisect.bisect_right(interval_ends, read_ref_start)
    last = bisect.bisect_left(interval_starts, read_ref_end)
    if last <= first:
        return range(0, 0)
    return range(first, last)


def _alignment_segment_key(read: pysam.AlignedSegment) -> tuple[Any, ...]:
    return (
        _read_key(read),
        int(read.reference_start),
        int(read.reference_end or read.reference_start),
        read.cigarstring or "",
        bool(read.is_supplementary),
        bool(read.is_reverse),
    )


def _record_from_projected(
    read: pysam.AlignedSegment,
    projected,
    window,
    read_sequence: str,
    aln_query_start: int,
) -> dict[str, Any] | None:
    query_start = int(projected.query_start) - int(aln_query_start)
    query_end = int(projected.query_end) - int(aln_query_start)
    if query_start < 0 or query_end <= query_start or query_end > len(read_sequence):
        return None
    read_start = int(projected.window_start) - int(window.raw_start)
    read_end = int(projected.window_end) - int(window.raw_start)
    if read_start < 0 or read_end <= read_start or read_end > int(window.actual_length):
        return None
    return {
        "read_key": _read_key(read),
        "shared_read_key": _alignment_segment_key(read),
        "read_sequence": read_sequence,
        "read_haplotype": _read_haplotype(read),
        "read_phase_set": _read_phase_set(read),
        "read_index": -1,
        "query_start": query_start,
        "query_end": query_end,
        "read_start": read_start,
        "read_end": read_end,
        "read_strand": -1 if read.is_reverse else 1,
        "reference_start": int(read.reference_start),
        "reference_end": int(read.reference_end or read.reference_start),
        "cigar": read.cigarstring or "",
        "cigar_ops_packed": _pack_cigartuples(read.cigartuples),
        "segment_type": 1 if read.is_supplementary else 0,
    }


def _project_identity_read_to_window_segment(
    read: pysam.AlignedSegment,
    *,
    window_start: int,
    window_end: int,
) -> ProjectedWindowSegment | None:
    """Fast path for SNP-only/identity HAP coordinates.

    In identity coordinate mode, individual raw coordinates are equal to BAM
    reference coordinates. We can therefore avoid dense coord.ref_to_raw lookup
    and project CIGAR operations directly against the target window.
    """
    if not read.query_sequence or not read.cigartuples:
        return None

    query_position = 0
    reference_position = int(read.reference_start)
    raw_min: int | None = None
    raw_max: int | None = None
    previous_raw: int | None = None
    selected_raw_min: int | None = None
    selected_raw_max: int | None = None
    selected_query_min: int | None = None
    selected_query_max: int | None = None

    def _select_query(q0: int, q1: int, raw0: int, raw1: int) -> None:
        nonlocal selected_query_min, selected_query_max
        nonlocal selected_raw_min, selected_raw_max
        if q1 <= q0:
            return
        selected_query_min = q0 if selected_query_min is None else min(selected_query_min, q0)
        selected_query_max = q1 - 1 if selected_query_max is None else max(selected_query_max, q1 - 1)
        selected_raw_min = raw0 if selected_raw_min is None else min(selected_raw_min, raw0)
        selected_raw_max = raw1 if selected_raw_max is None else max(selected_raw_max, raw1)

    for operation, length_value in read.cigartuples:
        length = int(length_value)
        if operation in (0, 7, 8):  # M, =, X
            block_ref_start = reference_position
            block_ref_end = reference_position + length
            if length > 0:
                raw_min = block_ref_start if raw_min is None else min(raw_min, block_ref_start)
                raw_max = block_ref_end - 1 if raw_max is None else max(raw_max, block_ref_end - 1)
                previous_raw = block_ref_end - 1
            overlap_start = max(block_ref_start, int(window_start))
            overlap_end = min(block_ref_end, int(window_end))
            if overlap_start < overlap_end:
                q0 = query_position + (overlap_start - block_ref_start)
                q1 = query_position + (overlap_end - block_ref_start)
                _select_query(q0, q1, overlap_start, overlap_end - 1)
            query_position += length
            reference_position += length
        elif operation == 1:  # insertion relative to reference
            anchor = previous_raw + 1 if previous_raw is not None else reference_position
            if window_start <= anchor < window_end:
                _select_query(query_position, query_position + length, anchor, anchor)
            query_position += length
        elif operation in (2, 3):  # deletion or skipped reference
            reference_position += length
        elif operation == 4:  # soft clip
            query_position += length
        elif operation in (5, 6):  # hard clip or padding
            continue
        else:
            raise ValueError(f"Unsupported CIGAR operation: {operation}")

    if (
        raw_min is None
        or raw_max is None
        or selected_query_min is None
        or selected_query_max is None
        or selected_raw_min is None
        or selected_raw_max is None
    ):
        return None

    return ProjectedWindowSegment(
        raw_start=raw_min,
        raw_end=raw_max + 1,
        query_start=selected_query_min,
        query_end=selected_query_max + 1,
        strand="-" if read.is_reverse else "+",
        window_start=selected_raw_min,
        window_end=selected_raw_max + 1,
    )


def _can_reuse_identity_projection(selected: dict[str, dict[str, Any]]) -> bool:
    if len(selected) < 2:
        return False
    signatures = []
    for data in selected.values():
        coord = data.get("coord")
        if getattr(coord, "source", "") != "identity":
            return False
        signatures.append(tuple((int(w.raw_start), int(w.actual_length)) for w in data["windows"]))
    return all(signature == signatures[0] for signature in signatures[1:])


def _append_projected_record(
    *,
    data: dict[str, Any],
    wi: int,
    record: dict[str, Any],
    timing: dict[str, Any] | None,
) -> bool:
    segment_key = (
        record["shared_read_key"],
        int(record["query_start"]),
        int(record["query_end"]),
        int(record["read_start"]),
        int(record["read_end"]),
    )
    step_started = time.perf_counter()
    if segment_key in data["seen"][wi]:
        _timing_add(timing, "chunk_seen_segment_dedup", time.perf_counter() - step_started)
        return False
    data["seen"][wi].add(segment_key)
    _timing_add(timing, "chunk_seen_segment_dedup", time.perf_counter() - step_started)
    data["records"][wi].append(record)
    data["counts"]["kept_segments"] += 1
    return True


def build_paired_hap_window_artifacts(
    *,
    sample: str,
    contig: str,
    bam_path: Path,
    hap_contexts: dict[str, dict[str, Any]],
    reference_fasta: Path | None = None,
    max_windows: int | None = None,
    window_start_index: int = 0,
    config: ExportConfig = ExportConfig(),
    timing: dict[str, Any] | None = None,
) -> dict[str, dict[str, Any]]:
    """Build HAP1/HAP2 chunk artifacts with one union BAM fetch pass.

    This is the optimized path for the shared-read schema. It selects the same
    eligible-window slice from every available HAP context, loads each HAP chain
    map for the chunk once, computes the union of all reference fetch intervals,
    fetches BAM alignments once, and projects the fetched alignments into each
    HAP/window coordinate system.
    """
    config.validate()
    if window_start_index < 0:
        raise ValueError("window_start_index must be non-negative")

    selected: dict[str, dict[str, Any]] = {}
    fetch_intervals: list[tuple[int, int]] = []
    for label, context in hap_contexts.items():
        step_started = time.perf_counter()
        eligible_rows = context["eligible_rows"]
        stop = None if max_windows is None else window_start_index + max_windows
        rows = eligible_rows[window_start_index:stop]
        if not rows:
            continue
        windows = [row[2] for row in rows]
        _timing_add(timing, "chunk_select_windows", time.perf_counter() - step_started)
        original = context["original"]
        flank = 1_000_000
        raw_range = (
            max(0, windows[0].raw_start - flank),
            min(len(original), windows[-1].raw_end + flank),
        )
        step_started = time.perf_counter()
        coord = _coordinate_map_from_context(context, raw_range)
        _timing_add(timing, "chunk_build_coordinate_maps", time.perf_counter() - step_started)
        step_started = time.perf_counter()
        window_ref_intervals = [reference_interval(coord, w.raw_start, w.raw_end) for w in windows]
        fetch_intervals.extend(window_ref_intervals)
        _timing_add(timing, "chunk_make_fetch_intervals", time.perf_counter() - step_started)
        selected[label] = {
            "rows": rows,
            "windows": windows,
            "coord": coord,
            "window_ref_intervals": window_ref_intervals,
            "window_ref_starts": [int(x[0]) for x in window_ref_intervals],
            "window_ref_ends": [int(x[1]) for x in window_ref_intervals],
            "records": [[] for _ in windows],
            "seen": [set() for _ in windows],
            "counts": Counter(),
        }
    if not selected:
        raise ValueError(
            f"no coordinate-eligible window selected: sample={sample} contig={contig} start={window_start_index}"
        )

    step_started = time.perf_counter()
    merged_fetch_intervals = _merge_intervals(fetch_intervals)
    _timing_add(timing, "chunk_merge_fetch_intervals", time.perf_counter() - step_started)
    _timing_set_count(timing, "fetch_interval_count", len(fetch_intervals))
    _timing_set_count(timing, "merged_fetch_interval_count", len(merged_fetch_intervals))

    fetched_reads: dict[tuple[Any, ...], pysam.AlignedSegment] = {}
    common_counts = Counter()
    raw_fetch_count = 0
    step_started = time.perf_counter()
    with _alignment_file(bam_path, reference_fasta) as bam:
        bam_contig = contig if contig in bam.references else contig.removeprefix("chr")
        for r0, r1 in merged_fetch_intervals:
            for read in bam.fetch(bam_contig, r0, r1):
                raw_fetch_count += 1
                reason = rejection_reason(read, config)
                if reason:
                    common_counts[f"excluded_{reason}"] += 1
                    continue
                key = _alignment_segment_key(read)
                fetched_reads.setdefault(key, read)
    _timing_add(timing, "chunk_bam_fetch_filter", time.perf_counter() - step_started)
    _timing_set_count(timing, "raw_fetch_read_count", raw_fetch_count)
    _timing_set_count(timing, "deduplicated_fetched_read_count", len(fetched_reads))

    identity_reuse = _can_reuse_identity_projection(selected)
    _timing_set_count(timing, "identity_reuse_enabled", 1 if identity_reuse else 0)
    projection_attempt_count = 0
    record_count = 0
    projection_reuse_record_count = 0
    labels = list(selected)
    reuse_base_label = labels[0] if identity_reuse else ""
    # 입력 순서를 보존하면서 처리한 BAM 객체와 디코딩 캐시를 즉시 해제함.
    # 투영 결과와 원본 객체를 청크 전체에 걸쳐 동시에 보관하지 않음.
    for read_key in tuple(fetched_reads):
        read = fetched_reads.pop(read_key)
        step_started = time.perf_counter()
        read_sequence, aln_query_start, _aln_query_end = _clipping_removed_query_sequence(read)
        _timing_add(timing, "chunk_prepare_read_sequence", time.perf_counter() - step_started)
        if not read_sequence:
            for data in selected.values():
                data["counts"]["excluded_empty_mapped_query"] += 1
            continue
        read_ref_start = int(read.reference_start)
        read_ref_end = int(read.reference_end or read.reference_start)

        if identity_reuse:
            base_data = selected[reuse_base_label]
            step_started = time.perf_counter()
            candidate_indices = _candidate_window_indices(
                interval_starts=base_data["window_ref_starts"],
                interval_ends=base_data["window_ref_ends"],
                read_ref_start=read_ref_start,
                read_ref_end=read_ref_end,
            )
            candidate_indices = list(candidate_indices)
            _timing_add(timing, "chunk_candidate_window_lookup", time.perf_counter() - step_started)
            projection_attempt_count += len(candidate_indices)
            projection_reuse_record_count += len(candidate_indices) * max(0, len(labels) - 1)
            for wi in candidate_indices:
                window = base_data["windows"][wi]
                step_started = time.perf_counter()
                projected = _project_identity_read_to_window_segment(
                    read,
                    window_start=window.raw_start,
                    window_end=window.raw_start + window.actual_length,
                )
                _timing_add(timing, "chunk_project_reads", time.perf_counter() - step_started)
                if projected is None:
                    for data in selected.values():
                        data["counts"]["excluded_unprojectable"] += 1
                    continue
                if projected.raw_start >= window.raw_end or projected.raw_end <= window.raw_start:
                    continue
                step_started = time.perf_counter()
                record = _record_from_projected(read, projected, window, read_sequence, aln_query_start)
                _timing_add(timing, "chunk_record_build", time.perf_counter() - step_started)
                if record is None:
                    for data in selected.values():
                        data["counts"]["excluded_invalid_segment_coordinates"] += 1
                    continue
                for label in labels:
                    data = selected[label]
                    target_window = data["windows"][wi]
                    target_record = record if label == reuse_base_label else dict(record)
                    if target_window.raw_start != window.raw_start or target_window.actual_length != window.actual_length:
                        target_record = dict(record)
                        target_record["read_start"] = int(projected.window_start) - int(target_window.raw_start)
                        target_record["read_end"] = int(projected.window_end) - int(target_window.raw_start)
                    if _append_projected_record(data=data, wi=wi, record=target_record, timing=timing):
                        record_count += 1
            continue

        for label, data in selected.items():
            coord = data["coord"]
            step_started = time.perf_counter()
            candidate_indices = _candidate_window_indices(
                interval_starts=data["window_ref_starts"],
                interval_ends=data["window_ref_ends"],
                read_ref_start=read_ref_start,
                read_ref_end=read_ref_end,
            )
            candidate_indices = list(candidate_indices)
            _timing_add(timing, "chunk_candidate_window_lookup", time.perf_counter() - step_started)
            projection_attempt_count += len(candidate_indices)
            for wi in candidate_indices:
                window = data["windows"][wi]
                step_started = time.perf_counter()
                if getattr(coord, "source", "") == "identity":
                    projected = _project_identity_read_to_window_segment(
                        read,
                        window_start=window.raw_start,
                        window_end=window.raw_start + window.actual_length,
                    )
                elif getattr(coord, "source", "") == "chain":
                    projected = project_read_to_window_segment_blockwise(
                        read,
                        coord,
                        window.raw_start,
                        window.raw_start + window.actual_length,
                    )
                else:
                    projected = project_read_to_window_segment(
                        read,
                        coord,
                        window.raw_start,
                        window.raw_start + window.actual_length,
                    )
                _timing_add(timing, "chunk_project_reads", time.perf_counter() - step_started)
                if projected is None:
                    data["counts"]["excluded_unprojectable"] += 1
                    continue
                if projected.raw_start >= window.raw_end or projected.raw_end <= window.raw_start:
                    continue
                step_started = time.perf_counter()
                record = _record_from_projected(read, projected, window, read_sequence, aln_query_start)
                _timing_add(timing, "chunk_record_build", time.perf_counter() - step_started)
                if record is None:
                    data["counts"]["excluded_invalid_segment_coordinates"] += 1
                    continue
                if _append_projected_record(data=data, wi=wi, record=record, timing=timing):
                    record_count += 1

    _timing_set_count(timing, "projection_attempt_count", projection_attempt_count)
    _timing_set_count(timing, "projection_reuse_record_count", projection_reuse_record_count)
    _timing_set_count(timing, "record_count", record_count)
    for key, value in common_counts.items():
        _timing_set_count(timing, str(key), int(value))

    result: dict[str, dict[str, Any]] = {}
    for label, data in selected.items():
        rows = data["rows"]
        context = hap_contexts[label]
        counts = Counter(common_counts)
        counts.update(data["counts"])
        attrs = {
            "sample": sample,
            "label": label,
            "contig": contig,
            "coordinate_system": "individual_raw_0based_halfopen",
            "terminal_n_left": context["summary"]["terminal_n_left"],
            "terminal_n_right": context["summary"]["terminal_n_right"],
            "random_offset": context["summary"]["random_offset"],
            "random_offset_unit": context["summary"]["random_offset_unit"],
            "valid_window_start_index": window_start_index,
            "sequence_valid_window_count_total": context["summary"]["sequence_valid_window_count_total"],
            "coordinate_eligible_window_count_total": context["summary"]["coordinate_eligible_window_count_total"],
            "coordinate_excluded_window_count": context["summary"]["coordinate_excluded_window_count"],
            "coordinate_excluded_bases": context["summary"]["coordinate_excluded_bases"],
            "selected_coordinate_mapped_bases": [row[0] for row in rows],
            "selected_coordinate_coverage": [row[1] for row in rows],
            "coordinate_map_source": context["summary"].get("coordinate_map_source", "chain"),
            "source_chain_path": str(context["chain_path"] or ""),
            "source_individual_fasta": str(context["individual_fasta"]),
            "source_bam_path": str(bam_path),
            "window_index_path": context["summary"].get("window_index_path", ""),
            "window_index_schema": context["summary"].get("window_index_schema", ""),
            "window_index_kept_count": context["summary"].get("window_index_kept_count", 0),
            "window_index_total_count": context["summary"].get("window_index_total_count", 0),
            "window_index_dropped_n": context["summary"].get("window_index_dropped_n", 0),
            "window_index_dropped_blacklist": context["summary"].get("window_index_dropped_blacklist", 0),
            "window_index_drop_reason_priority": context["summary"].get("window_index_drop_reason_priority", ""),
            "config": config.to_dict(),
        }
        result[label] = {
            "label": label,
            "windows": data["windows"],
            "records": data["records"],
            "counts": counts,
            "attrs": attrs,
            "window_count": len(data["windows"]),
            "segment_count": sum(len(x) for x in data["records"]),
        }
    return result


class SharedHapChrWriter:
    """Write one sample/contig zarr with shared reads and HAP-specific windows."""

    def __init__(
        self,
        path: Path,
        *,
        sample: str,
        contig: str,
        hap_window_counts: dict[str, int],
        window_size: int,
        attrs: dict | None = None,
    ) -> None:
        self.path = Path(path)
        self.sample = sample
        self.contig = contig
        self.window_size = int(window_size)
        self.hap_window_counts = {label: int(count) for label, count in hap_window_counts.items()}
        self.hap_window_cursor = {label: 0 for label in self.hap_window_counts}
        self.hap_segment_cursor = {label: 0 for label in self.hap_window_counts}
        self.hap_counts: dict[str, Counter[str]] = {label: Counter() for label in self.hap_window_counts}
        self.hap_coverage: dict[str, list[float]] = {label: [] for label in self.hap_window_counts}
        self.read_key_to_index: dict[tuple[str, int], int] = {}
        self.read_haplotype_by_index: list[int] = []
        self.read_phase_set_by_index: list[int] = []
        self.read_length_sum = 0
        self.read_length_min: int | None = None
        self.read_length_max: int | None = None
        self.haplotype_conflict_count = 0
        self.phase_set_conflict_count = 0

        comp = Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE)
        root = zarr.open_group(str(self.path), mode="w")
        root.attrs.update({
            "schema": SHARED_SCHEMA,
            "schema_version": SCHEMA_VERSION,
            "base_encoding": "literal_string",
            "sample": sample,
            "contig": contig,
            "coordinate_system": "haplotype_raw_0based_halfopen_per_HAP_group",
            "read_sequence_storage": READ_SEQUENCE_STORAGE,
            "read_haplotype_storage": READ_HAPLOTYPE_STORAGE,
            "read_phase_set_storage": READ_PHASE_SET_STORAGE,
            "read_cigar_storage": READ_CIGAR_STORAGE,
            "read_reference_storage": READ_REFERENCE_STORAGE,
            "phase_set_storage": READ_PHASE_SET_STORAGE,
            "window_size": self.window_size,
            "hap_labels": list(self.hap_window_counts),
            "created_at": utc_now(),
            **(attrs or {}),
        })
        self.root = root
        reads = root.create_group("reads")
        self.read_arrays = {
            "read_sequence": reads.empty("read_sequence", shape=(0,), dtype=object, object_codec=VLenUTF8(), chunks=(READ_CHUNK_SIZE,)),
            "read_sequence_offsets": reads.empty("read_sequence_offsets", shape=(1,), dtype=np.int64, chunks=(READ_CHUNK_SIZE,), compressor=comp),
            "read_length": reads.empty("read_length", shape=(0,), dtype=np.int32, chunks=(READ_CHUNK_SIZE,), compressor=comp),
            "read_haplotype": reads.empty("read_haplotype", shape=(0,), dtype=np.int8, chunks=(READ_CHUNK_SIZE,), compressor=comp),
            "read_phase_set": reads.empty("read_phase_set", shape=(0,), dtype=np.int64, chunks=(READ_CHUNK_SIZE,), compressor=comp),
            "read_reference_start": reads.empty("read_reference_start", shape=(0,), dtype=np.int64, chunks=(READ_CHUNK_SIZE,), compressor=comp),
            "read_reference_end": reads.empty("read_reference_end", shape=(0,), dtype=np.int64, chunks=(READ_CHUNK_SIZE,), compressor=comp),
            "read_segment_type": reads.empty("read_segment_type", shape=(0,), dtype=np.int8, chunks=(READ_CHUNK_SIZE,), compressor=comp),
            "read_cigar_offset": reads.empty("read_cigar_offset", shape=(0,), dtype=np.int64, chunks=(READ_CHUNK_SIZE,), compressor=comp),
            "read_cigar_len": reads.empty("read_cigar_len", shape=(0,), dtype=np.int32, chunks=(READ_CHUNK_SIZE,), compressor=comp),
        }
        self.read_arrays["read_sequence_offsets"][0] = 0
        cigar = root.create_group("cigar")
        self.cigar_ops = cigar.empty("ops_packed", shape=(0,), dtype=np.uint32, chunks=(SEGMENT_CHUNK_SIZE,), compressor=comp)
        self.cigar_op_count = 0
        self.hap_groups: dict[str, zarr.Group] = {}
        self.hap_arrays: dict[str, dict[str, Any]] = {}
        for label, total_windows in self.hap_window_counts.items():
            group = root.create_group(label)
            self.hap_groups[label] = group
            self.hap_arrays[label] = {
                "sequence": group.empty("sequence", shape=(total_windows,), dtype=object, object_codec=VLenUTF8(), chunks=(1,)),
                "raw_start": group.empty("raw_start", shape=(total_windows,), dtype=np.int64, chunks=(min(total_windows + 1, 1024),), compressor=comp),
                "actual_length": group.empty("actual_length", shape=(total_windows,), dtype=np.int64, chunks=(min(total_windows + 1, 1024),), compressor=comp),
                "read_index": group.empty("read_index", shape=(0,), dtype=np.int64, chunks=(SEGMENT_CHUNK_SIZE,), compressor=comp),
                "query_start": group.empty("query_start", shape=(0,), dtype=np.int32, chunks=(SEGMENT_CHUNK_SIZE,), compressor=comp),
                "query_end": group.empty("query_end", shape=(0,), dtype=np.int32, chunks=(SEGMENT_CHUNK_SIZE,), compressor=comp),
                "read_start": group.empty("read_start", shape=(0,), dtype=np.int64, chunks=(SEGMENT_CHUNK_SIZE,), compressor=comp),
                "read_end": group.empty("read_end", shape=(0,), dtype=np.int64, chunks=(SEGMENT_CHUNK_SIZE,), compressor=comp),
                "read_strand": group.empty("read_strand", shape=(0,), dtype=np.int8, chunks=(SEGMENT_CHUNK_SIZE,), compressor=comp),
                "window_read_offsets": group.empty("window_read_offsets", shape=(total_windows + 1,), dtype=np.int64, chunks=(min(total_windows + 1, 1024),), compressor=comp),
            }
            self.hap_arrays[label]["window_read_offsets"][0] = 0

    def _register_read(self, record: dict[str, Any]) -> int:
        # A shared read entry represents one BAM alignment segment, not just a
        # query name. Supplementary alignments can have different CIGAR/clipping
        # boundaries, so using read name alone would make query_start/query_end
        # ambiguous. The same segment is still shared across HAP1 and HAP2.
        key = tuple(record.get("shared_read_key", record["read_key"]))
        existing = self.read_key_to_index.get(key)
        hap = int(record["read_haplotype"])
        phase_set = int(record.get("read_phase_set", 0))
        if existing is not None:
            old_hap = self.read_haplotype_by_index[existing]
            if old_hap != hap and hap != 0:
                self.read_haplotype_by_index[existing] = 0
                self.read_arrays["read_haplotype"][existing] = 0
                self.haplotype_conflict_count += 1
            old_phase_set = self.read_phase_set_by_index[existing]
            if old_phase_set != phase_set and phase_set != 0:
                self.read_phase_set_by_index[existing] = 0
                self.read_arrays["read_phase_set"][existing] = 0
                self.phase_set_conflict_count += 1
            return existing

        index = len(self.read_key_to_index)
        sequence = str(record["read_sequence"]).upper()
        if set(sequence) - set("ACGTN"):
            raise ValueError("invalid read alphabet outside A/C/G/T/N")
        self.read_key_to_index[key] = index
        self.read_haplotype_by_index.append(hap)
        self.read_phase_set_by_index.append(phase_set)
        new = index + 1
        for name in ("read_sequence", "read_length", "read_haplotype", "read_phase_set", "read_reference_start", "read_reference_end", "read_segment_type", "read_cigar_offset", "read_cigar_len"):
            self.read_arrays[name].resize((new,))
        self.read_arrays["read_sequence_offsets"].resize((new + 1,))
        self.read_arrays["read_sequence"][index] = sequence
        self.read_arrays["read_length"][index] = len(sequence)
        self.read_arrays["read_haplotype"][index] = hap
        self.read_arrays["read_phase_set"][index] = phase_set
        cigar_ops = np.asarray(record.get("cigar_ops_packed", []), dtype=np.uint32)
        cigar_offset = self.cigar_op_count
        cigar_len = int(cigar_ops.shape[0])
        if cigar_len:
            self.cigar_ops.resize((cigar_offset + cigar_len,))
            self.cigar_ops[cigar_offset:cigar_offset + cigar_len] = cigar_ops
            self.cigar_op_count += cigar_len
        self.read_arrays["read_reference_start"][index] = int(record.get("reference_start", -1))
        self.read_arrays["read_reference_end"][index] = int(record.get("reference_end", -1))
        self.read_arrays["read_segment_type"][index] = int(record.get("segment_type", 0))
        self.read_arrays["read_cigar_offset"][index] = cigar_offset
        self.read_arrays["read_cigar_len"][index] = cigar_len
        previous_offset = int(self.read_arrays["read_sequence_offsets"][index])
        self.read_arrays["read_sequence_offsets"][new] = previous_offset + len(sequence)
        self.read_length_min = len(sequence) if self.read_length_min is None else min(self.read_length_min, len(sequence))
        self.read_length_max = len(sequence) if self.read_length_max is None else max(self.read_length_max, len(sequence))
        self.read_length_sum += len(sequence)
        return index

    def _register_reads_bulk(self, records: list[dict[str, Any]]) -> None:
        """Register all new shared alignment-segment sequences in one zarr write.

        Per-read zarr resize is very slow. This method first identifies all new
        shared_read_key values in a computation chunk, then appends their strings,
        lengths, haplotype labels, phase-set ids, and cumulative offsets in a
        single batch. Existing shared reads are left in place. If the same key is
        observed with conflicting non-zero HP or PS labels, the stored value is
        downgraded to 0 (unknown), which is safer than assigning it to a wrong
        haplotype or phase block.
        """
        pending: dict[tuple[Any, ...], dict[str, Any]] = {}
        pending_order: list[tuple[Any, ...]] = []
        for record in records:
            key = tuple(record.get("shared_read_key", record["read_key"]))
            hap = int(record["read_haplotype"])
            phase_set = int(record.get("read_phase_set", 0))
            existing = self.read_key_to_index.get(key)
            if existing is not None:
                old_hap = self.read_haplotype_by_index[existing]
                if old_hap != hap and hap != 0:
                    self.read_haplotype_by_index[existing] = 0
                    self.read_arrays["read_haplotype"][existing] = 0
                    self.haplotype_conflict_count += 1
                old_phase_set = self.read_phase_set_by_index[existing]
                if old_phase_set != phase_set and phase_set != 0:
                    self.read_phase_set_by_index[existing] = 0
                    self.read_arrays["read_phase_set"][existing] = 0
                    self.phase_set_conflict_count += 1
                continue
            if key in pending:
                old_hap = int(pending[key]["read_haplotype"])
                if old_hap != hap and hap != 0:
                    pending[key]["read_haplotype"] = 0
                old_phase_set = int(pending[key].get("read_phase_set", 0))
                if old_phase_set != phase_set and phase_set != 0:
                    pending[key]["read_phase_set"] = 0
                continue
            pending[key] = record
            pending_order.append(key)

        if not pending_order:
            return

        old = len(self.read_key_to_index)
        sequences = [str(pending[key]["read_sequence"]).upper() for key in pending_order]
        for sequence in sequences:
            if set(sequence) - set("ACGTN"):
                raise ValueError("invalid read alphabet outside A/C/G/T/N")
        lengths = np.asarray([len(sequence) for sequence in sequences], dtype=np.int32)
        haplotypes = np.asarray([int(pending[key]["read_haplotype"]) for key in pending_order], dtype=np.int8)
        phase_sets = np.asarray([int(pending[key].get("read_phase_set", 0)) for key in pending_order], dtype=np.int64)
        reference_starts = np.asarray([int(pending[key].get("reference_start", -1)) for key in pending_order], dtype=np.int64)
        reference_ends = np.asarray([int(pending[key].get("reference_end", -1)) for key in pending_order], dtype=np.int64)
        segment_types = np.asarray([int(pending[key].get("segment_type", 0)) for key in pending_order], dtype=np.int8)
        cigar_lists = [np.asarray(pending[key].get("cigar_ops_packed", []), dtype=np.uint32) for key in pending_order]
        cigar_lens = np.asarray([int(values.shape[0]) for values in cigar_lists], dtype=np.int32)
        cigar_offsets = self.cigar_op_count + np.concatenate(([0], np.cumsum(cigar_lens[:-1], dtype=np.int64))).astype(np.int64)
        if np.any(lengths <= 0):
            raise ValueError("empty shared read sequence")
        if np.any(phase_sets < 0):
            raise ValueError("negative phase-set id")
        if np.any(reference_ends < reference_starts):
            raise ValueError("invalid reference interval")
        if np.any(segment_types < 0):
            raise ValueError("negative segment type")
        new = old + len(pending_order)

        total_new_cigar_ops = int(cigar_lens.sum())
        if total_new_cigar_ops:
            old_cigar = self.cigar_op_count
            new_cigar = old_cigar + total_new_cigar_ops
            self.cigar_ops.resize((new_cigar,))
            self.cigar_ops[old_cigar:new_cigar] = np.concatenate(cigar_lists).astype(np.uint32, copy=False)
            self.cigar_op_count = new_cigar

        for name in ("read_sequence", "read_length", "read_haplotype", "read_phase_set", "read_reference_start", "read_reference_end", "read_segment_type", "read_cigar_offset", "read_cigar_len"):
            self.read_arrays[name].resize((new,))
        self.read_arrays["read_sequence_offsets"].resize((new + 1,))
        self.read_arrays["read_sequence"][old:new] = np.asarray(sequences, dtype=object)
        self.read_arrays["read_length"][old:new] = lengths
        self.read_arrays["read_haplotype"][old:new] = haplotypes
        self.read_arrays["read_phase_set"][old:new] = phase_sets
        self.read_arrays["read_reference_start"][old:new] = reference_starts
        self.read_arrays["read_reference_end"][old:new] = reference_ends
        self.read_arrays["read_segment_type"][old:new] = segment_types
        self.read_arrays["read_cigar_offset"][old:new] = cigar_offsets
        self.read_arrays["read_cigar_len"][old:new] = cigar_lens
        previous_offset = int(self.read_arrays["read_sequence_offsets"][old])
        self.read_arrays["read_sequence_offsets"][old + 1:new + 1] = previous_offset + np.cumsum(lengths, dtype=np.int64)

        for offset, key in enumerate(pending_order, start=old):
            self.read_key_to_index[key] = offset
            self.read_haplotype_by_index.append(int(haplotypes[offset - old]))
            self.read_phase_set_by_index.append(int(phase_sets[offset - old]))
        batch_min = int(lengths.min())
        batch_max = int(lengths.max())
        self.read_length_min = batch_min if self.read_length_min is None else min(self.read_length_min, batch_min)
        self.read_length_max = batch_max if self.read_length_max is None else max(self.read_length_max, batch_max)
        self.read_length_sum += int(lengths.sum())

    def append_hap_chunk(self, label: str, artifacts: dict) -> None:
        if label not in self.hap_arrays:
            raise ValueError(f"unknown hap label: {label}")
        windows = artifacts["windows"]
        records = artifacts["records"]
        if len(windows) != len(records):
            raise ValueError("windows/records length mismatch")
        if self.hap_window_cursor[label] + len(windows) > self.hap_window_counts[label]:
            raise ValueError(f"too many windows appended for {label}")

        flat_records = [record for window_records in records for record in window_records]
        self._register_reads_bulk(flat_records)
        self._write_hap_provenance_attrs(label, artifacts.get("attrs", {}))

        arrays = self.hap_arrays[label]
        wlo = self.hap_window_cursor[label]
        whi = wlo + len(windows)
        arrays["sequence"][wlo:whi] = np.asarray([w.sequence for w in windows], dtype=object)
        arrays["raw_start"][wlo:whi] = np.asarray([w.raw_start for w in windows], dtype=np.int64)
        arrays["actual_length"][wlo:whi] = np.asarray([w.actual_length for w in windows], dtype=np.int64)

        packed: list[tuple[int, int, int, int, int, int]] = []
        local_offsets = [0]
        for window, window_records in zip(windows, records):
            row_values: list[tuple[int, int, int, int, int, int]] = []
            for record in window_records:
                key = tuple(record.get("shared_read_key", record["read_key"]))
                read_index = self.read_key_to_index[key]
                query_start = int(record["query_start"])
                query_end = int(record["query_end"])
                read_start = int(record["read_start"])
                read_end = int(record["read_end"])
                read_strand = int(record["read_strand"])
                if query_start < 0 or query_end <= query_start:
                    raise ValueError("invalid query interval")
                if read_index < 0 or read_index >= len(self.read_key_to_index):
                    raise ValueError("shared read index out of range")
                read_length = int(self.read_arrays["read_length"][read_index])
                if query_end > read_length:
                    raise ValueError(f"query interval outside shared read: {query_start}-{query_end} length={read_length}")
                if read_start < 0 or read_end <= read_start or read_end > int(window.actual_length):
                    raise ValueError(f"invalid hap read interval: {read_start}-{read_end}")
                if read_strand not in (-1, 1):
                    raise ValueError("invalid read strand")
                row_values.append((read_index, query_start, query_end, read_start, read_end, read_strand))
            row_values.sort(key=lambda x: (x[0], x[1], x[2], x[3], x[4]))
            packed.extend(row_values)
            local_offsets.append(len(packed))

        old = self.hap_segment_cursor[label]
        new = old + len(packed)
        for name in ("read_index", "query_start", "query_end", "read_start", "read_end", "read_strand"):
            arrays[name].resize((new,))
        if packed:
            values = np.asarray(packed, dtype=np.int64)
            arrays["read_index"][old:new] = values[:, 0].astype(np.int64)
            arrays["query_start"][old:new] = values[:, 1].astype(np.int32)
            arrays["query_end"][old:new] = values[:, 2].astype(np.int32)
            arrays["read_start"][old:new] = values[:, 3].astype(np.int64)
            arrays["read_end"][old:new] = values[:, 4].astype(np.int64)
            arrays["read_strand"][old:new] = values[:, 5].astype(np.int8)
        arrays["window_read_offsets"][wlo + 1:whi + 1] = np.asarray(local_offsets[1:], dtype=np.int64) + old
        self.hap_window_cursor[label] = whi
        self.hap_segment_cursor[label] = new
        self.hap_counts[label].update(artifacts.get("counts", {}))
        self.hap_coverage[label].extend(float(x) for x in artifacts.get("attrs", {}).get("selected_coordinate_coverage", []))

    def _write_hap_provenance_attrs(self, label: str, attrs: dict[str, Any]) -> None:
        """Persist HAP-level coordinate provenance once per HAP group.

        Chunk-specific fields such as selected coordinate coverage are summarized
        at the root during finalize. The HAP group keeps only stable provenance
        needed after raw intermediates are compacted and removed.
        """
        group = self.hap_groups[label]
        for key in HAP_PROVENANCE_ATTR_KEYS:
            if key in attrs and key not in group.attrs:
                group.attrs[key] = attrs[key]

    def finalize(self, *, extra_attrs: dict | None = None, validation_level: str = "standard") -> dict:
        for label, total in self.hap_window_counts.items():
            if self.hap_window_cursor[label] != total:
                raise ValueError(f"window count mismatch for {label}: {self.hap_window_cursor[label]} != {total}")
        n_reads = len(self.read_key_to_index)
        self.root.attrs.update({
            "shared_read_count": n_reads,
            "shared_read_length_min": self.read_length_min,
            "shared_read_length_max": self.read_length_max,
            "shared_read_length_mean": (self.read_length_sum / n_reads) if n_reads else None,
            "read_haplotype_conflict_count": self.haplotype_conflict_count,
            "read_phase_set_conflict_count": self.phase_set_conflict_count,
            "shared_read_cigar_op_count": int(self.cigar_op_count),
            "read_cigar_storage": READ_CIGAR_STORAGE,
            "read_reference_storage": READ_REFERENCE_STORAGE,
            "hap_segment_count": {label: self.hap_segment_cursor[label] for label in self.hap_window_counts},
            "hap_filter_counts": {label: dict(sorted(self.hap_counts[label].items())) for label in self.hap_window_counts},
            "hap_coordinate_coverage_min": {label: min(self.hap_coverage[label]) if self.hap_coverage[label] else None for label in self.hap_window_counts},
            "hap_coordinate_coverage_max": {label: max(self.hap_coverage[label]) if self.hap_coverage[label] else None for label in self.hap_window_counts},
            "completed_at": utc_now(),
            **(extra_attrs or {}),
        })
        qc = validate_shared_hap_zarr(self.path, window_size=self.window_size, validation_level=validation_level)
        write_json(self.path / "chr_qc.json", qc)
        return qc


def validate_shared_hap_zarr(path: Path, *, window_size: int, validation_level: str = "standard") -> dict:
    """Validate a shared-HAP zarr.

    validation_level controls only expensive string scans; it does not change the
    preprocessing output.

    - full: check every read string and every window sequence/padding.
    - standard: skip expensive object-string scans, but keep all shape/range/
      offset/HP/PS/segment coordinate checks. This is the operational default.
    - fast: same string-scan policy as standard, intended for quick resume checks
      where a cheaper validation call is preferred over full content validation.
    """
    if validation_level not in VALIDATION_LEVELS:
        raise ValueError(f"validation_level must be one of: {', '.join(VALIDATION_LEVELS)}")
    root = zarr.open_group(str(path), mode="r")
    errors: list[str] = []
    schema = root.attrs.get("schema")
    if schema not in COMPATIBLE_SHARED_SCHEMAS:
        errors.append("unexpected_schema")
    if "reads" not in root:
        errors.append("missing_reads_group")
        return {"status": "fail", "errors": errors}
    reads = root["reads"]
    requires_cigar = root.attrs.get("schema_version") == SCHEMA_VERSION or bool(root.attrs.get("read_cigar_storage"))
    read_required = {"read_sequence", "read_sequence_offsets", "read_length", "read_haplotype", "read_phase_set"}
    if requires_cigar:
        read_required.update({"read_reference_start", "read_reference_end", "read_segment_type", "read_cigar_offset", "read_cigar_len"})
    if not read_required.issubset(set(reads.array_keys())):
        errors.append("missing_read_arrays")
    n_reads = int(reads["read_sequence"].shape[0]) if "read_sequence" in reads else 0
    if len(reads["read_sequence_offsets"]) != n_reads + 1:
        errors.append("invalid_read_sequence_offsets_length")
    read_lengths = reads["read_length"][:] if "read_length" in reads else np.asarray([], dtype=np.int32)
    if len(read_lengths) != n_reads:
        errors.append("read_length_mismatch")
    phase_sets = reads["read_phase_set"][:] if "read_phase_set" in reads else np.asarray([], dtype=np.int64)
    if len(phase_sets) != n_reads:
        errors.append("read_phase_set_length_mismatch")
    reference_starts = reads["read_reference_start"][:] if "read_reference_start" in reads else np.asarray([], dtype=np.int64)
    reference_ends = reads["read_reference_end"][:] if "read_reference_end" in reads else np.asarray([], dtype=np.int64)
    cigar_offsets = reads["read_cigar_offset"][:] if "read_cigar_offset" in reads else np.asarray([], dtype=np.int64)
    cigar_lens = reads["read_cigar_len"][:] if "read_cigar_len" in reads else np.asarray([], dtype=np.int32)
    cigar_ops = root["cigar"]["ops_packed"] if "cigar" in root and "ops_packed" in root["cigar"] else None
    if requires_cigar and (len(reference_starts) != n_reads or len(reference_ends) != n_reads):
        errors.append("read_reference_interval_length_mismatch")
    if requires_cigar and (len(cigar_offsets) != n_reads or len(cigar_lens) != n_reads or cigar_ops is None):
        errors.append("read_cigar_length_mismatch")
    if n_reads:
        if not set(np.unique(reads["read_haplotype"][:])).issubset({0, 1, 2}):
            errors.append("invalid_read_haplotype")
        if len(phase_sets) == n_reads and np.any(phase_sets < 0):
            errors.append("invalid_read_phase_set")
        if np.any(read_lengths <= 0):
            errors.append("invalid_read_length")
        if requires_cigar and len(reference_starts) == n_reads and len(reference_ends) == n_reads and np.any(reference_ends < reference_starts):
            errors.append("invalid_read_reference_interval")
        if requires_cigar and len(cigar_offsets) == n_reads and len(cigar_lens) == n_reads and cigar_ops is not None:
            if np.any(cigar_offsets < 0) or np.any(cigar_lens < 0) or np.any(cigar_offsets + cigar_lens > int(cigar_ops.shape[0])):
                errors.append("invalid_read_cigar_offsets")
    offsets = reads["read_sequence_offsets"][:] if "read_sequence_offsets" in reads else np.asarray([], dtype=np.int64)
    if len(offsets) and (int(offsets[0]) != 0 or np.any(np.diff(offsets) < 0)):
        errors.append("invalid_read_sequence_offsets_value")
    length_errors = 0
    alphabet_errors = 0
    if validation_level == "full":
        for lo in range(0, n_reads, READ_CHUNK_SIZE):
            hi = min(n_reads, lo + READ_CHUNK_SIZE)
            values = reads["read_sequence"][lo:hi]
            for j, value in enumerate(values, start=lo):
                seq = str(value)
                if len(seq) != int(read_lengths[j]):
                    length_errors += 1
                if set(seq) - set("ACGTN"):
                    alphabet_errors += 1
        if length_errors:
            errors.append("read_sequence_length_mismatch")
        if alphabet_errors:
            errors.append("invalid_read_sequence_alphabet")

    hap_qc: dict[str, dict[str, Any]] = {}
    for label in root.attrs.get("hap_labels", HAP_LABELS):
        if label not in root:
            errors.append(f"missing_{label}_group")
            continue
        group = root[label]
        required = {"sequence", "raw_start", "actual_length", "read_index", "query_start", "query_end", "read_start", "read_end", "read_strand", "window_read_offsets"}
        if not required.issubset(set(group.array_keys())):
            errors.append(f"missing_{label}_arrays")
            continue
        n_windows = int(group["sequence"].shape[0])
        n_segments = int(group["read_index"].shape[0])
        hap_errors: list[str] = []
        for name in ("query_start", "query_end", "read_start", "read_end", "read_strand"):
            if int(group[name].shape[0]) != n_segments:
                hap_errors.append(f"{name}_length_mismatch")
        window_offsets = group["window_read_offsets"][:]
        if len(window_offsets) != n_windows + 1:
            hap_errors.append("invalid_window_read_offsets_length")
        elif int(window_offsets[0]) != 0 or int(window_offsets[-1]) != n_segments or np.any(np.diff(window_offsets) < 0):
            hap_errors.append("invalid_window_read_offsets_values")
        actual = group["actual_length"][:]
        if n_windows and np.any(actual <= 0):
            hap_errors.append("invalid_actual_length")
        if n_segments:
            read_index = group["read_index"][:]
            query_start = group["query_start"][:]
            query_end = group["query_end"][:]
            read_start = group["read_start"][:]
            read_end = group["read_end"][:]
            read_strand = group["read_strand"][:]
            if int(read_index.min()) < 0 or int(read_index.max()) >= n_reads:
                hap_errors.append("read_index_out_of_range")
            else:
                if np.any(query_start < 0) or np.any(query_start >= query_end) or np.any(query_end > read_lengths[read_index]):
                    hap_errors.append("query_interval_out_of_range")
            if np.any(read_start < 0) or np.any(read_start >= read_end):
                hap_errors.append("invalid_read_interval")
            if not set(np.unique(read_strand)).issubset({-1, 1}):
                hap_errors.append("invalid_read_strand")
        for wi in range(n_windows):
            actual_length = int(actual[wi])
            if validation_level == "full":
                seq = str(group["sequence"][wi])
                if len(seq) != window_size or set(seq[:actual_length]) - set("ACGTN") or seq[actual_length:] != "P" * (window_size - actual_length):
                    hap_errors.append("invalid_window_sequence_or_padding")
                    break
            lo, hi = int(window_offsets[wi]), int(window_offsets[wi + 1])
            if hi > lo and np.any(group["read_end"][lo:hi] > actual_length):
                hap_errors.append("read_interval_outside_window")
                break
        coverage = [float(x) for x in root.attrs.get("hap_coordinate_coverage", {}).get(label, [])]
        hap_qc[label] = {
            "errors": hap_errors,
            "window_count": n_windows,
            "segment_count": n_segments,
            "coordinate_coverage": {
                "min": min(coverage) if coverage else root.attrs.get("hap_coordinate_coverage_min", {}).get(label),
                "max": max(coverage) if coverage else root.attrs.get("hap_coordinate_coverage_max", {}).get(label),
            },
        }
        errors.extend(f"{label}:{err}" for err in hap_errors)

    return {
        "status": "pass" if not errors else "fail",
        "errors": errors,
        "sample": root.attrs.get("sample"),
        "contig": root.attrs.get("contig"),
        "schema": root.attrs.get("schema"),
        "validation_level": validation_level,
        "shared_read_count": n_reads,
        "shared_read_length": {
            "min": int(np.min(read_lengths)) if n_reads else None,
            "max": int(np.max(read_lengths)) if n_reads else None,
            "mean": float(np.mean(read_lengths)) if n_reads else None,
        },
        "read_cigar": {
            "required": bool(requires_cigar),
            "op_count": int(cigar_ops.shape[0]) if cigar_ops is not None else 0,
            "storage": root.attrs.get("read_cigar_storage"),
        },
        "read_haplotype_counts": {
            str(int(value)): int(count)
            for value, count in zip(*np.unique(reads["read_haplotype"][:] if n_reads else np.asarray([], dtype=np.int8), return_counts=True))
        },
        "read_phase_set": {
            "present_count": int(np.count_nonzero(phase_sets)) if len(phase_sets) == n_reads else 0,
            "missing_count": int(n_reads - np.count_nonzero(phase_sets)) if len(phase_sets) == n_reads else n_reads,
            "unique_nonzero_count": int(len(np.unique(phase_sets[phase_sets != 0]))) if len(phase_sets) == n_reads else 0,
            "conflict_count": int(root.attrs.get("read_phase_set_conflict_count", 0)),
        },
        "haplotypes": hap_qc,
    }
