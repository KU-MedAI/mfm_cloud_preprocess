from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pysam
import zarr

from .config import ExportConfig
from .coordinates import (
    chain_query_mapped_intervals,
    interval_mapping_coverage,
    load_chain_map,
    load_paf_map,
    paf_query_mapped_intervals,
    project_read_to_window,
    reference_interval,
)
from .filters import rejection_reason
from .minimap2_backend import run_minimap2
from .storage import overlapping_alignment_indices, write_store
from .model_storage import write_model_store
from .windowing import deterministic_offset, make_windows, trim_terminal_n


def _fasta_sequence(path: Path) -> tuple[str, str]:
    with pysam.FastaFile(str(path)) as fasta:
        if len(fasta.references) != 1:
            raise ValueError("individual FASTA must contain exactly one contig")
        name = fasta.references[0]
        return name, fasta.fetch(name).upper()


def summarize_windows(*, sample: str, contig: str, individual_fasta: Path,
                      chain_path: Path,
                      config: ExportConfig = ExportConfig()) -> dict:
    config.validate()
    query_name, original = _fasta_sequence(individual_fasta)
    trimmed = trim_terminal_n(original)
    offset_unit = config.effective_offset_unit()
    offset = deterministic_offset(config.seed, sample, contig, offset_unit)
    windows = make_windows(
        trimmed.sequence,
        source_start=trimmed.source_start,
        offset=offset,
        window_size=config.window_size,
        max_n_ratio=config.max_n_ratio,
    )
    mapped_intervals = chain_query_mapped_intervals(
        chain_path, query_name=query_name, target_name=contig
    )
    rows = [
        (*interval_mapping_coverage(mapped_intervals, window.raw_start, window.raw_end), window)
        for window in windows
    ]
    eligible = [row for row in rows if row[1] >= config.min_coordinate_coverage]
    coverages = [float(row[1]) for row in rows]
    return {
        "sample": sample,
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
        "coordinate_coverage_min": min(coverages) if coverages else None,
        "coordinate_coverage_max": max(coverages) if coverages else None,
        "coordinate_coverage_threshold": config.min_coordinate_coverage,
    }


def build_readset_artifacts(*, sample: str, contig: str, individual_fasta: Path,
                            bam_path: Path,
                            chain_path: Path | None = None,
                            coordinate_backend: str = "chain",
                            reference_fasta: Path | None = None,
                            paf_path: Path | None = None,
                            minimap2: str = "minimap2",
                            max_windows: int | None = None,
                            window_start_index: int = 0,
                            config: ExportConfig = ExportConfig()) -> dict:
    """Build one chunk of model-input records without writing shard files."""
    config.validate()
    if window_start_index < 0:
        raise ValueError("window_start_index must be non-negative")
    query_name, original = _fasta_sequence(individual_fasta)
    trimmed = trim_terminal_n(original)
    ref_name = contig
    ref_len = None
    if coordinate_backend == "chain":
        if chain_path is None or not chain_path.exists() or not chain_path.stat().st_size:
            raise FileNotFoundError(f"chain file is required: {chain_path}")
    elif coordinate_backend == "minimap2-paf":
        if reference_fasta is None:
            raise ValueError("reference_fasta is required for minimap2-paf backend")
        if paf_path is None:
            raise ValueError("paf_path is required for minimap2-paf backend")
        with pysam.FastaFile(str(reference_fasta)) as ref:
            ref_name = contig if contig in ref.references else contig.removeprefix("chr")
            if ref_name not in ref.references:
                raise ValueError(f"{contig} absent from reference FASTA")
            ref_len = ref.get_reference_length(ref_name)
        if not paf_path.exists() or not paf_path.stat().st_size:
            run_minimap2(reference_fasta, individual_fasta, paf_path, minimap2)
    else:
        raise ValueError(f"unsupported coordinate_backend: {coordinate_backend}")

    offset_unit = config.effective_offset_unit()
    offset = deterministic_offset(config.seed, sample, contig, offset_unit)
    all_windows = make_windows(
        trimmed.sequence,
        source_start=trimmed.source_start,
        offset=offset,
        window_size=config.window_size,
        max_n_ratio=config.max_n_ratio,
    )
    if coordinate_backend == "chain":
        mapped_intervals = chain_query_mapped_intervals(
            chain_path, query_name=query_name, target_name=contig
        )
    else:
        mapped_intervals = paf_query_mapped_intervals(
            paf_path, query_name=query_name, target_name=ref_name
        )
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
            f"no coordinate-eligible window selected: start={window_start_index}, "
            f"sequence_valid={len(all_windows)}, eligible={len(eligible_rows)}, "
            f"threshold={config.min_coordinate_coverage}"
        )

    flank = 1_000_000
    raw_range = (
        max(0, windows[0].raw_start - flank),
        min(len(original), windows[-1].raw_end + flank),
    )
    if coordinate_backend == "chain":
        coord = load_chain_map(
            chain_path,
            len(original),
            query_name=query_name,
            target_name=contig,
            raw_range=raw_range,
        )
    else:
        coord = load_paf_map(
            paf_path,
            len(original),
            int(ref_len),
            query_name=query_name,
            target_name=ref_name,
            raw_range=raw_range,
        )

    reads: list[str] = []
    alignments: list[dict] = []
    associations: list[list[int]] = [[] for _ in windows]
    association_sets: list[set[int]] = [set() for _ in windows]
    model_records: list[list[dict]] = [[] for _ in windows]
    counts = Counter()
    name_to_read: dict[tuple[str, int], int] = {}
    segment_to_alignment: dict[tuple, int] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        bam_contig = contig if contig in bam.references else contig.removeprefix("chr")
        for wi, window in enumerate(windows):
            r0, r1 = reference_interval(coord, window.raw_start, window.raw_end)
            for read in bam.fetch(bam_contig, r0, r1):
                reason = rejection_reason(read, config)
                if reason:
                    counts[f"excluded_{reason}"] += 1
                    continue
                projected = project_read_to_window(
                    read, coord, window.raw_start, window.raw_start + window.actual_length
                )
                if projected is None:
                    counts["excluded_unprojectable"] += 1
                    continue
                raw_start, raw_end = projected.raw_start, projected.raw_end
                query_start, query_end = projected.query_start, projected.query_end
                strand = projected.strand
                if raw_start >= window.raw_end or raw_end <= window.raw_start:
                    continue
                if (
                    not projected.sequence
                    or projected.window_start is None
                    or projected.window_end is None
                ):
                    counts["excluded_no_query_base_in_window"] += 1
                    continue
                mate_code = 1 if read.is_read1 else (2 if read.is_read2 else 0)
                read_key = (read.query_name, mate_code)
                key = (
                    read_key,
                    read.reference_start,
                    read.cigarstring,
                    read.is_supplementary,
                    read.is_reverse,
                )
                if key in segment_to_alignment:
                    ai = segment_to_alignment[key]
                else:
                    ri = name_to_read.get(read_key)
                    if ri is None:
                        ri = len(reads)
                        name_to_read[read_key] = ri
                        reads.append(read.query_sequence.upper())
                    ai = len(alignments)
                    segment_to_alignment[key] = ai
                    alignments.append({
                        "key": key,
                        "read_index": ri,
                        "raw_start": raw_start,
                        "raw_end": raw_end,
                        "query_start": query_start,
                        "query_end": query_end,
                        "reference_start": int(read.reference_start),
                        "reference_end": int(read.reference_end or read.reference_start),
                        "strand": -1 if strand == "-" else 1,
                        "segment_type": 1 if read.is_supplementary else 0,
                        "mapq": int(read.mapping_quality),
                        "duplicate": bool(read.is_duplicate),
                        "cigar": read.cigarstring or "",
                    })
                    counts["kept_segments"] += 1
                if ai is not None and ai not in association_sets[wi]:
                    association_sets[wi].add(ai)
                    associations[wi].append(ai)
                    model_records[wi].append({
                        "sequence": projected.sequence,
                        "start": projected.window_start - window.raw_start,
                        "end": projected.window_end - window.raw_start,
                        "strand": int(alignments[ai]["strand"]),
                        "alignment_index": int(ai),
                    })

    attrs = {
        "sample": sample,
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
        "coordinate_map_source": coordinate_backend,
        "source_chain_path": str(chain_path) if chain_path else None,
        "source_paf_path": str(paf_path) if paf_path else None,
        "source_reference_fasta": str(reference_fasta) if reference_fasta else None,
        "config": config.to_dict(),
        "source_bam_path": str(bam_path),
        "source_bam_size": bam_path.stat().st_size,
        "source_bam_sha256_prefix": _hash_prefix(bam_path),
    }
    return {
        "windows": windows,
        "reads": reads,
        "alignments": alignments,
        "associations": associations,
        "model_records": model_records,
        "counts": counts,
        "attrs": attrs,
        "window_count": len(windows),
        "read_count": len(reads),
        "alignment_count": len(alignments),
        "model_read_count": sum(len(records) for records in model_records),
    }


def export_readset(*, sample: str, contig: str, individual_fasta: Path,
                   bam_path: Path, output_dir: Path,
                   chain_path: Path | None = None,
                   coordinate_backend: str = "chain",
                   reference_fasta: Path | None = None,
                   paf_path: Path | None = None,
                   minimap2: str = "minimap2",
                   max_windows: int | None = None,
                   window_start_index: int = 0,
                   config: ExportConfig = ExportConfig()) -> dict:
    config.validate()
    if window_start_index < 0:
        raise ValueError("window_start_index must be non-negative")
    output_dir.mkdir(parents=True, exist_ok=True)
    query_name, original = _fasta_sequence(individual_fasta)
    trimmed = trim_terminal_n(original)
    ref_name = contig
    ref_len = None
    if coordinate_backend == "chain":
        if chain_path is None or not chain_path.exists() or not chain_path.stat().st_size:
            raise FileNotFoundError(f"chain file is required: {chain_path}")
    elif coordinate_backend == "minimap2-paf":
        if reference_fasta is None:
            raise ValueError("reference_fasta is required for minimap2-paf backend")
        if paf_path is None:
            raise ValueError("paf_path is required for minimap2-paf backend")
        with pysam.FastaFile(str(reference_fasta)) as ref:
            ref_name = contig if contig in ref.references else contig.removeprefix("chr")
            if ref_name not in ref.references:
                raise ValueError(f"{contig} absent from reference FASTA")
            ref_len = ref.get_reference_length(ref_name)
        if not paf_path.exists() or not paf_path.stat().st_size:
            run_minimap2(reference_fasta, individual_fasta, paf_path, minimap2)
    else:
        raise ValueError(f"unsupported coordinate_backend: {coordinate_backend}")
    offset_unit = config.effective_offset_unit()
    offset = deterministic_offset(config.seed, sample, contig, offset_unit)
    all_windows = make_windows(trimmed.sequence, source_start=trimmed.source_start,
                               offset=offset, window_size=config.window_size,
                               max_n_ratio=config.max_n_ratio)
    if coordinate_backend == "chain":
        mapped_intervals = chain_query_mapped_intervals(
            chain_path, query_name=query_name, target_name=contig
        )
    else:
        mapped_intervals = paf_query_mapped_intervals(
            paf_path, query_name=query_name, target_name=ref_name
        )
    coverage_rows = [
        (*interval_mapping_coverage(mapped_intervals, window.raw_start, window.raw_end), window)
        for window in all_windows
    ]
    eligible_rows = [
        row for row in coverage_rows if row[1] >= config.min_coordinate_coverage
    ]
    stop = None if max_windows is None else window_start_index + max_windows
    selected_rows = eligible_rows[window_start_index:stop]
    windows = [row[2] for row in selected_rows]
    if not windows:
        raise ValueError(
            f"no coordinate-eligible window selected: start={window_start_index}, "
            f"sequence_valid={len(all_windows)}, eligible={len(eligible_rows)}, "
            f"threshold={config.min_coordinate_coverage}"
        )
    # Include flanking individual coordinates so reads crossing a window boundary
    # retain their full projected interval before model-level overlap clipping.
    flank = 1_000_000
    raw_range = (max(0, windows[0].raw_start - flank),
                 min(len(original), windows[-1].raw_end + flank))
    if coordinate_backend == "chain":
        coord = load_chain_map(chain_path, len(original), query_name=query_name,
                               target_name=contig, raw_range=raw_range)
    else:
        coord = load_paf_map(paf_path, len(original), int(ref_len),
                             query_name=query_name, target_name=ref_name,
                             raw_range=raw_range)
    reads: list[str] = []
    alignments: list[dict] = []
    associations: list[list[int]] = [[] for _ in windows]
    association_sets: list[set[int]] = [set() for _ in windows]
    model_records: list[list[dict]] = [[] for _ in windows]
    counts = Counter()
    name_to_read: dict[tuple[str, int], int] = {}
    segment_to_alignment: dict[tuple, int] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        bam_contig = contig if contig in bam.references else contig.removeprefix("chr")
        for wi, window in enumerate(windows):
            r0, r1 = reference_interval(coord, window.raw_start, window.raw_end)
            for read in bam.fetch(bam_contig, r0, r1):
                reason = rejection_reason(read, config)
                if reason:
                    counts[f"excluded_{reason}"] += 1
                    continue
                projected = project_read_to_window(
                    read, coord, window.raw_start, window.raw_start + window.actual_length
                )
                if projected is None:
                    counts["excluded_unprojectable"] += 1
                    continue
                raw_start, raw_end = projected.raw_start, projected.raw_end
                query_start, query_end = projected.query_start, projected.query_end
                strand = projected.strand
                if raw_start >= window.raw_end or raw_end <= window.raw_start:
                    continue
                if (
                    not projected.sequence
                    or projected.window_start is None
                    or projected.window_end is None
                ):
                    counts["excluded_no_query_base_in_window"] += 1
                    continue
                mate_code = 1 if read.is_read1 else (2 if read.is_read2 else 0)
                read_key = (read.query_name, mate_code)
                key = (read_key, read.reference_start, read.cigarstring,
                       read.is_supplementary, read.is_reverse)
                if key in segment_to_alignment:
                    ai = segment_to_alignment[key]
                else:
                    ri = name_to_read.get(read_key)
                    if ri is None:
                        ri = len(reads); name_to_read[read_key] = ri
                        reads.append(read.query_sequence.upper())
                    ai = len(alignments)
                    segment_to_alignment[key] = ai
                    alignments.append({"key": key, "read_index": ri, "raw_start": raw_start,
                        "raw_end": raw_end, "query_start": query_start, "query_end": query_end,
                        "reference_start": int(read.reference_start),
                        "reference_end": int(read.reference_end or read.reference_start),
                        "strand": -1 if strand == "-" else 1,
                        "segment_type": 1 if read.is_supplementary else 0,
                        "mapq": int(read.mapping_quality), "duplicate": bool(read.is_duplicate),
                        "cigar": read.cigarstring or ""})
                    counts["kept_segments"] += 1
                if ai is not None and ai not in association_sets[wi]:
                    association_sets[wi].add(ai)
                    associations[wi].append(ai)
                    model_records[wi].append({
                        "sequence": projected.sequence,
                        "start": projected.window_start - window.raw_start,
                        "end": projected.window_end - window.raw_start,
                        "strand": int(alignments[ai]["strand"]),
                        "alignment_index": int(ai),
                    })
    attrs = {"sample": sample, "contig": contig, "coordinate_system": "individual_raw_0based_halfopen",
             "terminal_n_left": trimmed.source_start, "terminal_n_right": len(original)-trimmed.source_end,
             "random_offset": offset, "random_offset_unit": offset_unit,
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
             "coordinate_map_source": coordinate_backend,
             "source_chain_path": str(chain_path) if chain_path else None,
             "source_paf_path": str(paf_path) if paf_path else None,
             "source_reference_fasta": str(reference_fasta) if reference_fasta else None,
             "config": config.to_dict(),
             "source_bam_path": str(bam_path),
             "source_bam_size": bam_path.stat().st_size,
             "source_bam_sha256_prefix": _hash_prefix(bam_path)}
    store_path = output_dir / "readset.zarr"
    write_store(store_path, windows, reads, alignments, associations, attrs)
    write_model_store(output_dir / "model_input.zarr", windows, reads, alignments, associations,
                      {"sample": sample, "contig": contig, "random_offset": offset,
                       "random_offset_unit": offset_unit,
                       "coordinate_map_source": coordinate_backend,
                       "valid_window_start_index": window_start_index,
                       "sequence_valid_window_count_total": len(all_windows),
                       "coordinate_eligible_window_count_total": len(eligible_rows),
                       "coordinate_coverage_threshold": config.min_coordinate_coverage,
                       "selected_coordinate_coverage": [row[1] for row in selected_rows]},
                      model_records=model_records)
    qc = validate_and_report(store_path, output_dir / "model_input.zarr",
                             alignments, counts, output_dir, config.window_size)
    return qc


def _hash_prefix(path: Path, size: int = 1 << 20) -> str:
    with path.open("rb") as handle:
        return hashlib.sha256(handle.read(size)).hexdigest()


def validate_and_report(
    store_path: Path,
    model_store_path: Path,
    audit_rows: list[dict],
    counts: Counter,
    output: Path,
    window_size: int,
    audit_max_rows_per_window: int = 20_000,
) -> dict:
    root = zarr.open_group(str(store_path), mode="r")
    model = zarr.open_group(str(model_store_path), mode="r")
    errors = []
    n_windows = len(root["windows/raw_start"])
    n_reads = len(root["reads/sequence"])
    n_aln = len(root["alignments/read_index"])
    n_model = len(model["read_sequence"])
    coordinate_coverage = list(root.attrs.get("selected_coordinate_coverage", []))
    coordinate_threshold = float(
        root.attrs.get("config", {}).get("min_coordinate_coverage", 0.0)
    )
    if len(coordinate_coverage) != n_windows:
        errors.append("coordinate_coverage_count_mismatch")
    elif any(float(value) < coordinate_threshold for value in coordinate_coverage):
        errors.append("coordinate_coverage_below_threshold")

    if n_aln and int(root["alignments/read_index"][:].max()) >= n_reads:
        errors.append("read_index_out_of_range")
    if n_aln and np.any(root["alignments/raw_start"][:] >= root["alignments/raw_end"][:]):
        errors.append("invalid_alignment_interval")
    if n_aln:
        query_start = root["alignments/query_start"][:]
        query_end = root["alignments/query_end"][:]
        read_lengths = np.asarray([len(x) for x in root["reads/sequence"][:]], dtype=np.int64)
        read_index = root["alignments/read_index"][:]
        if np.any(query_start < 0) or np.any(query_start >= query_end) or np.any(query_end > read_lengths[read_index]):
            errors.append("invalid_query_interval")

    raw_offsets = root["window_alignments/offsets"][:]
    raw_indices = root["window_alignments/alignment_indices"][:]
    if (
        len(raw_offsets) != n_windows + 1
        or int(raw_offsets[0]) != 0
        or int(raw_offsets[-1]) != len(raw_indices)
        or np.any(np.diff(raw_offsets) < 0)
    ):
        errors.append("invalid_raw_window_offsets")
    if len(raw_indices) and (int(raw_indices.min()) < 0 or int(raw_indices.max()) >= n_aln):
        errors.append("raw_alignment_index_out_of_range")

    required_model_arrays = {
        "sequence", "raw_start", "actual_length", "read_sequence",
        "read_sequence_offsets", "read_start", "read_end", "read_length",
        "read_strand", "alignment_index", "window_read_offsets",
    }
    if not required_model_arrays.issubset(set(model.array_keys())):
        errors.append("missing_model_array")
    if len(model["sequence"]) != n_windows:
        errors.append("model_window_count_mismatch")
    if not np.array_equal(model["raw_start"][:], root["windows/raw_start"][:]):
        errors.append("model_window_start_mismatch")
    if not np.array_equal(model["actual_length"][:], root["windows/actual_length"][:]):
        errors.append("model_actual_length_mismatch")

    model_offsets = model["window_read_offsets"][:]
    seq_offsets = model["read_sequence_offsets"][:]
    model_fields = ["read_start", "read_end", "read_length", "read_strand", "alignment_index"]
    if any(len(model[name]) != n_model for name in model_fields):
        errors.append("model_record_field_length_mismatch")
    if (
        len(model_offsets) != n_windows + 1
        or int(model_offsets[0]) != 0
        or int(model_offsets[-1]) != n_model
        or np.any(np.diff(model_offsets) < 0)
    ):
        errors.append("invalid_model_window_offsets")
    if (
        len(seq_offsets) != n_model + 1
        or not np.array_equal(seq_offsets, np.arange(n_model + 1, dtype=np.int64))
    ):
        errors.append("invalid_read_sequence_offsets")

    model_alignment_indices = model["alignment_index"][:]
    if len(model_alignment_indices) and (
        int(model_alignment_indices.min()) < 0 or int(model_alignment_indices.max()) >= n_aln
    ):
        errors.append("model_alignment_index_out_of_range")
    starts = model["read_start"][:]
    ends = model["read_end"][:]
    lengths = model["read_length"][:]
    strands = model["read_strand"][:]
    if n_model and (np.any(starts < 0) or np.any(starts >= ends)):
        errors.append("invalid_model_read_interval")
    if n_model and not set(np.unique(strands)).issubset({-1, 1}):
        errors.append("invalid_model_strand")

    model_read_sequences = [str(value) for value in model["read_sequence"][:]]
    model_sequence_length_errors = 0
    model_alphabet_errors = 0
    for index, sequence in enumerate(model_read_sequences):
        if len(sequence) != int(lengths[index]):
            model_sequence_length_errors += 1
        if set(sequence) - set("ACGTN"):
            model_alphabet_errors += 1
    if model_sequence_length_errors:
        errors.append("model_read_length_mismatch")
    if model_alphabet_errors:
        errors.append("invalid_model_read_alphabet")

    window_sequence_errors = 0
    model_bounds_errors = 0
    association_errors = 0
    for wi in range(n_windows):
        sequence = str(model["sequence"][wi])
        actual = int(model["actual_length"][wi])
        if (
            len(sequence) != window_size
            or set(sequence[:actual]) - set("ACGTN")
            or not any(base in "ACGT" for base in sequence[:actual])
            or sequence[actual:] != "P" * (window_size - actual)
        ):
            window_sequence_errors += 1
        m0, m1 = int(model_offsets[wi]), int(model_offsets[wi + 1])
        if m1 > m0 and np.any(ends[m0:m1] > actual):
            model_bounds_errors += 1
        r0, r1 = int(raw_offsets[wi]), int(raw_offsets[wi + 1])
        if not np.array_equal(model_alignment_indices[m0:m1], raw_indices[r0:r1]):
            association_errors += 1
    if window_sequence_errors:
        errors.append("invalid_window_sequence_or_padding")
    if model_bounds_errors:
        errors.append("model_read_outside_actual_window")
    if association_errors:
        errors.append("raw_model_association_mismatch")

    crop_counts = [len(overlapping_alignment_indices(root, i, 0)) for i in range(n_windows)]
    read_sequences = [str(value) for value in root["reads/sequence"][:]]
    columns = ["window_index", "window_raw_start", "window_raw_end", "alignment_index",
               "read_index", "raw_start", "raw_end", "query_start", "query_end", "strand",
               "segment_type", "read_length", "mapq_qc", "duplicate_qc", "bam_cigar_qc", "sequence_preview"]
    audit_rows_written = 0
    audit_sampled_windows = 0
    with (output / "alignment_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns); writer.writeheader()
        for wi in range(n_windows):
            indices = raw_indices[raw_offsets[wi]:raw_offsets[wi+1]]
            sorted_indices = sorted(indices, key=lambda i: (audit_rows[int(i)]["raw_start"], audit_rows[int(i)]["raw_end"]))
            if len(sorted_indices) > audit_max_rows_per_window:
                audit_sampled_windows += 1
                positions = np.linspace(
                    0, len(sorted_indices) - 1, audit_max_rows_per_window, dtype=np.int64
                )
                sorted_indices = [sorted_indices[int(pos)] for pos in positions]
            for ai in sorted_indices:
                a = audit_rows[int(ai)]; seq = read_sequences[a["read_index"]]
                writer.writerow({"window_index": wi, "window_raw_start": int(root["windows/raw_start"][wi]),
                    "window_raw_end": int(root["windows/raw_start"][wi])+int(root["windows/actual_length"][wi]),
                    "alignment_index": int(ai), "read_index": a["read_index"], "raw_start": a["raw_start"],
                    "raw_end": a["raw_end"], "query_start": a["query_start"], "query_end": a["query_end"],
                    "strand": "+" if a["strand"] == 1 else "-", "segment_type": "supplementary" if a["segment_type"] else "primary",
                    "read_length": len(seq), "mapq_qc": a["mapq"], "duplicate_qc": a["duplicate"],
                    "bam_cigar_qc": a["cigar"], "sequence_preview": seq[:40]})
                audit_rows_written += 1

    raw_read_lengths = np.asarray([len(x) for x in read_sequences], dtype=np.int64)
    qc = {
        "status": "pass" if not errors else "fail",
        "errors": errors,
        "window_count": n_windows,
        "read_count": n_reads,
        "alignment_count": n_aln,
        "model_read_count": n_model,
        "first_8k_alignment_counts": crop_counts,
        "raw_read_length": {
            "min": int(raw_read_lengths.min()) if n_reads else 0,
            "median": float(np.median(raw_read_lengths)) if n_reads else 0.0,
            "max": int(raw_read_lengths.max()) if n_reads else 0,
        },
        "model_read_length": {
            "min": int(lengths.min()) if n_model else 0,
            "median": float(np.median(lengths)) if n_model else 0.0,
            "max": int(lengths.max()) if n_model else 0,
        },
        "audit": {
            "rows_written": audit_rows_written,
            "max_rows_per_window": audit_max_rows_per_window,
            "sampled_window_count": audit_sampled_windows,
            "is_sampled": bool(audit_sampled_windows),
        },
        "coordinate_mapping": {
            "threshold": coordinate_threshold,
            "selected_coverage": [float(value) for value in coordinate_coverage],
            "sequence_valid_window_count_total": int(
                root.attrs.get("sequence_valid_window_count_total", n_windows)
            ),
            "eligible_window_count_total": int(
                root.attrs.get("coordinate_eligible_window_count_total", n_windows)
            ),
            "excluded_window_count": int(
                root.attrs.get("coordinate_excluded_window_count", 0)
            ),
            "excluded_bases": int(root.attrs.get("coordinate_excluded_bases", 0)),
        },
        "filter_counts": dict(sorted(counts.items())),
    }
    (output / "qc.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")
    return qc
