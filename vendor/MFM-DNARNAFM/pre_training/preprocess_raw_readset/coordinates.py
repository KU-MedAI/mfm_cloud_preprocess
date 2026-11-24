from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
import numpy as np

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

@dataclass(frozen=True)
class CoordinateMap:
    raw_to_ref: np.ndarray
    ref_to_raw: np.ndarray
    raw_offset: int
    ref_offset: int
    source: str = "unknown"
    ref_block_starts: np.ndarray | None = None
    ref_block_ends: np.ndarray | None = None
    ref_block_raw_starts: np.ndarray | None = None


@dataclass(frozen=True)
class ChainAlignment:
    """Parsed plus/plus UCSC chain alignment reused across window chunks.

    Chain parsing is pure text processing and can become visible overhead when a
    chromosome is split into many computation chunks. The preprocessing worker
    parses each HAP chain once in `prepare_hap_context()`, then builds small
    raw-range lookup arrays from these cached blocks for each chunk.
    """

    header: dict
    blocks: tuple[tuple[int, int, int], ...]


@dataclass(frozen=True)
class ProjectedWindowRead:
    raw_start: int
    raw_end: int
    query_start: int
    query_end: int
    strand: str
    sequence: str
    window_start: int | None
    window_end: int | None


@dataclass(frozen=True)
class ProjectedWindowSegment:
    """Window-local read segment coordinates without storing the segment bases.

    `query_start/query_end` are 0-based half-open coordinates on the original BAM
    query sequence. The caller may subtract `read.query_alignment_start` when the
    shared read storage keeps clipping-removed mapped query sequence only.

    `window_start/window_end` are individual-FASTA raw coordinates for the part
    of the read that overlaps the target window. They are not yet shifted by the
    window start.
    """

    raw_start: int
    raw_end: int
    query_start: int
    query_end: int
    strand: str
    window_start: int
    window_end: int

def _tag(fields: list[str], prefix: str) -> str | None:
    return next((x[len(prefix):] for x in fields[12:] if x.startswith(prefix)), None)


def paf_query_mapped_intervals(
    path: Path,
    *,
    query_name: str | None = None,
    target_name: str | None = None,
) -> list[tuple[int, int]]:
    """Return a merged union of plus-strand primary query intervals with M/=/X."""
    intervals: list[tuple[int, int]] = []
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if (
                len(fields) < 12
                or (query_name and fields[0] != query_name)
                or (target_name and fields[5] != target_name)
                or fields[4] != "+"
                or (_tag(fields, "tp:A:") not in (None, "P"))
            ):
                continue
            cigar = _tag(fields, "cg:Z:")
            if not cigar:
                raise ValueError("PAF must contain cg:Z CIGAR; run minimap2 with -c")
            query_position = int(fields[2])
            for length_text, operation in _CIGAR.findall(cigar):
                length = int(length_text)
                if operation in "M=X":
                    intervals.append((query_position, query_position + length))
                    query_position += length
                elif operation in "IS":
                    query_position += length
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        previous_start, previous_end = merged[-1]
        if start <= previous_end:
            merged[-1] = (previous_start, max(previous_end, end))
        else:
            merged.append((start, end))
    return merged


def interval_mapping_coverage(
    intervals: list[tuple[int, int]],
    start: int,
    end: int,
) -> tuple[int, float]:
    """Return mapped query bases and fraction for one half-open query interval."""
    if end <= start:
        return 0, 0.0
    mapped = 0
    for interval_start, interval_end in intervals:
        if interval_end <= start:
            continue
        if interval_start >= end:
            break
        mapped += max(0, min(end, interval_end) - max(start, interval_start))
    return mapped, mapped / (end - start)


def _chain_blocks(
    path: Path,
    *,
    query_name: str | None = None,
    target_name: str | None = None,
) -> tuple[dict, list[tuple[int, int, int]]]:
    """Parse a plus-strand UCSC chain as target(reference)->query(individual)."""
    header: dict | None = None
    blocks: list[tuple[int, int, int]] = []
    with Path(path).open() as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            fields = text.split()
            if fields[0] == "chain":
                if header is not None:
                    raise ValueError(f"multiple chains are not supported: {path}")
                if len(fields) != 13:
                    raise ValueError(f"invalid chain header: {line.rstrip()}")
                header = {
                    "score": int(fields[1]),
                    "target_name": fields[2],
                    "target_size": int(fields[3]),
                    "target_strand": fields[4],
                    "target_start": int(fields[5]),
                    "target_end": int(fields[6]),
                    "query_name": fields[7],
                    "query_size": int(fields[8]),
                    "query_strand": fields[9],
                    "query_start": int(fields[10]),
                    "query_end": int(fields[11]),
                    "chain_id": fields[12],
                }
                if target_name and header["target_name"] != target_name:
                    raise ValueError(
                        f"chain target mismatch: expected={target_name}, found={header['target_name']}"
                    )
                if query_name and header["query_name"] != query_name:
                    raise ValueError(
                        f"chain query mismatch: expected={query_name}, found={header['query_name']}"
                    )
                if header["target_strand"] != "+" or header["query_strand"] != "+":
                    raise ValueError("only plus/plus bcftools consensus chains are supported")
                t = header["target_start"]
                q = header["query_start"]
                continue
            if header is None:
                raise ValueError(f"chain data appeared before header: {path}")
            values = [int(value) for value in fields]
            if len(values) not in (1, 3):
                raise ValueError(f"invalid chain block: {line.rstrip()}")
            size = values[0]
            if size:
                blocks.append((q, t, size))
            q += size
            t += size
            if len(values) == 3:
                t += values[1]
                q += values[2]
    if header is None:
        raise ValueError(f"missing chain header: {path}")
    if not blocks:
        raise ValueError(f"no chain alignment blocks: {path}")
    return header, blocks


def load_chain_alignment(
    path: Path,
    *,
    query_name: str | None = None,
    target_name: str | None = None,
) -> ChainAlignment:
    """Parse a chain file once and return reusable alignment blocks."""
    header, blocks = _chain_blocks(path, query_name=query_name, target_name=target_name)
    return ChainAlignment(header=header, blocks=tuple(blocks))


def chain_query_mapped_intervals(
    path: Path,
    *,
    query_name: str | None = None,
    target_name: str | None = None,
) -> list[tuple[int, int]]:
    _header, blocks = _chain_blocks(
        path, query_name=query_name, target_name=target_name
    )
    intervals = sorted((q, q + length) for q, _t, length in blocks)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        previous_start, previous_end = merged[-1]
        if start <= previous_end:
            merged[-1] = (previous_start, max(previous_end, end))
        else:
            merged.append((start, end))
    return merged


def load_chain_map(
    path: Path,
    raw_length: int,
    *,
    query_name: str | None = None,
    target_name: str | None = None,
    raw_range: tuple[int, int] | None = None,
) -> CoordinateMap:
    """Load a bcftools consensus chain into individual<->reference lookup arrays."""
    alignment = load_chain_alignment(path, query_name=query_name, target_name=target_name)
    return chain_alignment_to_map(alignment, raw_length, raw_range=raw_range)


def identity_coordinate_map(
    raw_length: int,
    *,
    raw_range: tuple[int, int] | None = None,
) -> CoordinateMap:
    """Return an identity individual<->reference coordinate map.

    This is the correct coordinate model for SNP-only personal FASTA. SNP-only
    consensus changes bases but never changes sequence length, so individual
    coordinates and reference coordinates remain identical. No chain file should
    be created or required for that policy.
    """
    raw_lo, raw_hi = raw_range or (0, raw_length)
    if raw_lo < 0 or raw_hi < raw_lo or raw_hi > raw_length:
        raise ValueError(f"invalid identity coordinate range: {raw_lo}-{raw_hi} / {raw_length}")
    values = np.arange(raw_lo, raw_hi, dtype=np.int64)
    return CoordinateMap(
        raw_to_ref=values.copy(),
        ref_to_raw=values.copy(),
        raw_offset=raw_lo,
        ref_offset=raw_lo,
        source="identity",
        ref_block_starts=np.asarray([raw_lo], dtype=np.int64),
        ref_block_ends=np.asarray([raw_hi], dtype=np.int64),
        ref_block_raw_starts=np.asarray([raw_lo], dtype=np.int64),
    )


def chain_alignment_to_map(
    alignment: ChainAlignment,
    raw_length: int,
    *,
    raw_range: tuple[int, int] | None = None,
) -> CoordinateMap:
    """Build individual<->reference lookup arrays from cached chain blocks."""
    header = alignment.header
    if raw_length != header["query_size"]:
        raise ValueError(
            f"chain query length mismatch: fasta={raw_length}, chain={header['query_size']}"
        )
    raw_lo, raw_hi = raw_range or (0, raw_length)
    blocks = []
    for q, t, length in alignment.blocks:
        lo, hi = max(q, raw_lo), min(q + length, raw_hi)
        if lo < hi:
            blocks.append((lo, t + lo - q, hi - lo))
    if not blocks:
        raise ValueError("no chain mapping in requested raw range")
    ref_lo = min(t for _q, t, _n in blocks)
    ref_hi = max(t + n for _q, t, n in blocks)
    raw_to_ref = np.full(raw_hi - raw_lo, -1, dtype=np.int64)
    ref_to_raw = np.full(ref_hi - ref_lo, -1, dtype=np.int64)
    ref_block_starts: list[int] = []
    ref_block_ends: list[int] = []
    ref_block_raw_starts: list[int] = []
    for q, t, length in blocks:
        raw_to_ref[q - raw_lo:q - raw_lo + length] = np.arange(
            t, t + length, dtype=np.int64
        )
        ref_to_raw[t - ref_lo:t - ref_lo + length] = np.arange(
            q, q + length, dtype=np.int64
        )
        ref_block_starts.append(int(t))
        ref_block_ends.append(int(t + length))
        ref_block_raw_starts.append(int(q))
    return CoordinateMap(
        raw_to_ref,
        ref_to_raw,
        raw_lo,
        ref_lo,
        source="chain",
        ref_block_starts=np.asarray(ref_block_starts, dtype=np.int64),
        ref_block_ends=np.asarray(ref_block_ends, dtype=np.int64),
        ref_block_raw_starts=np.asarray(ref_block_raw_starts, dtype=np.int64),
    )

def load_paf_map(path: Path, raw_length: int, reference_length: int,
                 query_name: str | None = None, target_name: str | None = None,
                 raw_range: tuple[int, int] | None = None) -> CoordinateMap:
    raw_lo, raw_hi = raw_range or (0, raw_length)
    blocks: list[tuple[int, int, int]] = []
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) < 12 or (query_name and f[0] != query_name) or (target_name and f[5] != target_name):
                continue
            if f[4] != "+" or (_tag(f, "tp:A:") not in (None, "P")):
                continue
            cigar = _tag(f, "cg:Z:")
            if not cigar:
                raise ValueError("PAF must contain cg:Z CIGAR; run minimap2 with -c")
            q, r = int(f[2]), int(f[7])
            for length_s, op in _CIGAR.findall(cigar):
                length = int(length_s)
                if op in "M=X":
                    lo, hi = max(q, raw_lo), min(q + length, raw_hi)
                    if lo < hi:
                        blocks.append((lo, r + lo - q, hi - lo))
                    q += length; r += length
                elif op in "IS": q += length
                elif op in "DN": r += length
    if not blocks:
        raise ValueError(f"no primary plus-strand coordinate mapping in {path}")
    ref_lo = min(r for _, r, _ in blocks)
    ref_hi = max(r + n for _, r, n in blocks)
    raw_to_ref = np.full(raw_hi - raw_lo, -1, dtype=np.int32)
    ref_to_raw = np.full(ref_hi - ref_lo, -1, dtype=np.int32)
    ref_block_starts: list[int] = []
    ref_block_ends: list[int] = []
    ref_block_raw_starts: list[int] = []
    for q, r, length in blocks:
        raw_to_ref[q-raw_lo:q-raw_lo+length] = np.arange(r, r+length, dtype=np.int32)
        ref_to_raw[r-ref_lo:r-ref_lo+length] = np.arange(q, q+length, dtype=np.int32)
        ref_block_starts.append(int(r))
        ref_block_ends.append(int(r + length))
        ref_block_raw_starts.append(int(q))
    return CoordinateMap(
        raw_to_ref,
        ref_to_raw,
        raw_lo,
        ref_lo,
        source="paf",
        ref_block_starts=np.asarray(ref_block_starts, dtype=np.int64),
        ref_block_ends=np.asarray(ref_block_ends, dtype=np.int64),
        ref_block_raw_starts=np.asarray(ref_block_raw_starts, dtype=np.int64),
    )

def reference_interval(coord: CoordinateMap, raw_start: int, raw_end: int) -> tuple[int, int]:
    values = coord.raw_to_ref[raw_start-coord.raw_offset:raw_end-coord.raw_offset]
    values = values[values >= 0]
    if not len(values): raise ValueError("window has no reference mapping")
    return int(values.min()), int(values.max()) + 1

def project_read(read, coord: CoordinateMap) -> tuple[int, int, int, int, str] | None:
    mapped = []
    for q, r in read.get_aligned_pairs(matches_only=False):
        local = -1 if r is None else r - coord.ref_offset
        if q is not None and 0 <= local < len(coord.ref_to_raw) and coord.ref_to_raw[local] >= 0:
            mapped.append((q, int(coord.ref_to_raw[local])))
    if not mapped: return None
    query = [x[0] for x in mapped]; raw = [x[1] for x in mapped]
    return min(raw), max(raw)+1, min(query), max(query)+1, "-" if read.is_reverse else "+"

def project_read_pairs(read, coord: CoordinateMap) -> list[tuple[int, int]]:
    """Return (query_index, individual_raw_pos) pairs for aligned query bases.

    Soft-clipped query bases are absent from aligned_pairs. Query insertions have no
    reference position; they are anchored immediately after the previous mapped base
    so their sequence is retained in the window-local model representation.
    """
    pairs = read.get_aligned_pairs(matches_only=False)
    mapped: list[tuple[int, int | None]] = []
    for q, r in pairs:
        if q is None:
            continue
        raw = None
        if r is not None:
            local = int(r) - coord.ref_offset
            if 0 <= local < len(coord.ref_to_raw) and coord.ref_to_raw[local] >= 0:
                raw = int(coord.ref_to_raw[local])
        mapped.append((int(q), raw))
    result: list[tuple[int, int]] = []
    previous = None
    for i, (q, raw) in enumerate(mapped):
        if raw is None:
            following = next((x for _, x in mapped[i + 1:] if x is not None), None)
            raw = previous + 1 if previous is not None else (following if following is not None else None)
        if raw is not None:
            result.append((q, int(raw)))
            previous = int(raw)
    return result


def _sequence_runs(sequence: str, query_start: int, mask: np.ndarray) -> list[str]:
    """Return sequence slices for contiguous true runs without per-base strings."""
    selected = np.flatnonzero(mask)
    if not len(selected):
        return []
    breaks = np.flatnonzero(np.diff(selected) != 1) + 1
    return [
        sequence[query_start + int(run[0]):query_start + int(run[-1]) + 1]
        for run in np.split(selected, breaks)
    ]


def project_read_to_window(
    read,
    coord: CoordinateMap,
    window_start: int,
    window_end: int,
) -> ProjectedWindowRead | None:
    """Project one BAM segment and clip its query sequence to one individual window.

    CIGAR is traversed once. M/=/X and insertion query bases are retained,
    deletion/skipped-reference operations advance only the reference, and
    soft/hard clips are explicitly excluded. This avoids materializing every
    aligned pair and prevents soft clips from being mistaken for insertions.
    """
    if not read.query_sequence or not read.cigartuples:
        return None
    sequence = read.query_sequence.upper()
    query_position = 0
    reference_position = int(read.reference_start)
    raw_min: int | None = None
    raw_max: int | None = None
    previous_raw: int | None = None
    selected_min: int | None = None
    selected_max: int | None = None
    fragments: list[str] = []

    for operation, length in read.cigartuples:
        length = int(length)
        if operation in (0, 7, 8):  # M, =, X
            local = reference_position - coord.ref_offset
            block_start = max(0, -local)
            block_end = min(length, len(coord.ref_to_raw) - local)
            if block_start < block_end:
                values = coord.ref_to_raw[
                    local + block_start:local + block_end
                ]
                valid = values >= 0
                if np.any(valid):
                    valid_values = values[valid]
                    block_min = int(valid_values.min())
                    block_max = int(valid_values.max())
                    raw_min = block_min if raw_min is None else min(raw_min, block_min)
                    raw_max = block_max if raw_max is None else max(raw_max, block_max)
                    previous_raw = int(valid_values[-1])
                    selected = valid & (values >= window_start) & (values < window_end)
                    if np.any(selected):
                        selected_values = values[selected]
                        block_selected_min = int(selected_values.min())
                        block_selected_max = int(selected_values.max())
                        selected_min = (
                            block_selected_min
                            if selected_min is None else min(selected_min, block_selected_min)
                        )
                        selected_max = (
                            block_selected_max
                            if selected_max is None else max(selected_max, block_selected_max)
                        )
                        fragments.extend(_sequence_runs(
                            sequence,
                            query_position + block_start,
                            selected,
                        ))
            query_position += length
            reference_position += length
        elif operation == 1:  # insertion relative to GRCh38
            anchor = None
            if previous_raw is not None:
                anchor = previous_raw + 1
            else:
                local = reference_position - coord.ref_offset
                if 0 <= local < len(coord.ref_to_raw) and coord.ref_to_raw[local] >= 0:
                    anchor = int(coord.ref_to_raw[local])
            if anchor is not None and window_start <= anchor < window_end:
                fragments.append(sequence[query_position:query_position + length])
                selected_min = anchor if selected_min is None else min(selected_min, anchor)
                selected_max = anchor if selected_max is None else max(selected_max, anchor)
            query_position += length
        elif operation in (2, 3):  # deletion, skipped reference
            reference_position += length
        elif operation == 4:  # soft clip: query exists but is not mapped
            query_position += length
        elif operation in (5, 6):  # hard clip, padding
            continue
        else:
            raise ValueError(f"Unsupported CIGAR operation: {operation}")

    if raw_min is None or raw_max is None:
        return None
    return ProjectedWindowRead(
        raw_start=raw_min,
        raw_end=raw_max + 1,
        query_start=int(read.query_alignment_start or 0),
        query_end=int(read.query_alignment_end or 0),
        strand="-" if read.is_reverse else "+",
        sequence="".join(fragments),
        window_start=selected_min,
        window_end=None if selected_max is None else selected_max + 1,
    )


def project_read_to_window_segment(
    read,
    coord: CoordinateMap,
    window_start: int,
    window_end: int,
) -> ProjectedWindowSegment | None:
    """Project one BAM segment into a haplotype window as query coordinates.

    This is the shared-read-storage companion to `project_read_to_window`.
    It performs the expensive CIGAR/reference->individual coordinate traversal
    during preprocessing, but it does not materialize or store the clipped read
    segment sequence. Instead it returns the query interval that a dataloader can
    slice from the shared clipping-removed read sequence.

    The returned query interval is the smallest contiguous query range covering
    query bases selected for the window. M/=/X bases are selected when their
    projected individual coordinate lies inside the window. Insertions are kept
    when their anchor lies inside the window. Deletions/skipped reference
    positions do not contribute query bases. Soft/hard clips and padding are not
    selected.
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

    for operation, length in read.cigartuples:
        length = int(length)
        if operation in (0, 7, 8):  # M, =, X
            local = reference_position - coord.ref_offset
            block_start = max(0, -local)
            block_end = min(length, len(coord.ref_to_raw) - local)
            if block_start < block_end:
                values = coord.ref_to_raw[local + block_start:local + block_end]
                valid = values >= 0
                if np.any(valid):
                    valid_values = values[valid]
                    block_min = int(valid_values.min())
                    block_max = int(valid_values.max())
                    raw_min = block_min if raw_min is None else min(raw_min, block_min)
                    raw_max = block_max if raw_max is None else max(raw_max, block_max)
                    previous_raw = int(valid_values[-1])

                    selected = valid & (values >= window_start) & (values < window_end)
                    if np.any(selected):
                        selected_indices = np.flatnonzero(selected)
                        breaks = np.flatnonzero(np.diff(selected_indices) != 1) + 1
                        for run in np.split(selected_indices, breaks):
                            q0 = query_position + block_start + int(run[0])
                            q1 = query_position + block_start + int(run[-1]) + 1
                            selected_values = values[run]
                            _select_query(q0, q1, int(selected_values.min()), int(selected_values.max()))
            query_position += length
            reference_position += length
        elif operation == 1:  # insertion relative to GRCh38
            anchor = None
            if previous_raw is not None:
                anchor = previous_raw + 1
            else:
                local = reference_position - coord.ref_offset
                if 0 <= local < len(coord.ref_to_raw) and coord.ref_to_raw[local] >= 0:
                    anchor = int(coord.ref_to_raw[local])
            if anchor is not None and window_start <= anchor < window_end:
                _select_query(query_position, query_position + length, anchor, anchor)
            query_position += length
        elif operation in (2, 3):  # deletion, skipped reference
            reference_position += length
        elif operation == 4:  # soft clip
            query_position += length
        elif operation in (5, 6):  # hard clip, padding
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



def _iter_ref_block_overlaps(
    coord: CoordinateMap,
    ref_start: int,
    ref_end: int,
):
    starts = coord.ref_block_starts
    ends = coord.ref_block_ends
    raw_starts = coord.ref_block_raw_starts
    if starts is None or ends is None or raw_starts is None or not len(starts):
        return
    index = int(np.searchsorted(ends, int(ref_start), side="right"))
    while index < len(starts) and int(starts[index]) < int(ref_end):
        block_ref_start = int(starts[index])
        block_ref_end = int(ends[index])
        overlap_ref_start = max(int(ref_start), block_ref_start)
        overlap_ref_end = min(int(ref_end), block_ref_end)
        if overlap_ref_start < overlap_ref_end:
            block_raw_start = int(raw_starts[index])
            overlap_raw_start = block_raw_start + (overlap_ref_start - block_ref_start)
            yield overlap_ref_start, overlap_ref_end, overlap_raw_start
        index += 1


def _map_ref_position_to_raw_blockwise(coord: CoordinateMap, ref_position: int) -> int | None:
    starts = coord.ref_block_starts
    ends = coord.ref_block_ends
    raw_starts = coord.ref_block_raw_starts
    if starts is None or ends is None or raw_starts is None or not len(starts):
        return None
    index = int(np.searchsorted(ends, int(ref_position), side="right"))
    if index >= len(starts):
        return None
    block_ref_start = int(starts[index])
    block_ref_end = int(ends[index])
    if block_ref_start <= int(ref_position) < block_ref_end:
        return int(raw_starts[index]) + (int(ref_position) - block_ref_start)
    return None


def project_read_to_window_segment_blockwise(
    read,
    coord: CoordinateMap,
    window_start: int,
    window_end: int,
) -> ProjectedWindowSegment | None:
    # Chain fast path for WGS-only SNV+indel coordinates. It intersects CIGAR
    # M/=/X blocks with plus/plus chain blocks and avoids dense per-base masks.
    if coord.ref_block_starts is None or coord.ref_block_ends is None or coord.ref_block_raw_starts is None:
        return project_read_to_window_segment(read, coord, window_start, window_end)
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
            for overlap_ref_start, overlap_ref_end, overlap_raw_start in _iter_ref_block_overlaps(coord, block_ref_start, block_ref_end):
                overlap_len = overlap_ref_end - overlap_ref_start
                overlap_raw_end = overlap_raw_start + overlap_len
                raw_min = overlap_raw_start if raw_min is None else min(raw_min, overlap_raw_start)
                raw_max = overlap_raw_end - 1 if raw_max is None else max(raw_max, overlap_raw_end - 1)
                previous_raw = overlap_raw_end - 1

                selected_raw_start = max(overlap_raw_start, int(window_start))
                selected_raw_end = min(overlap_raw_end, int(window_end))
                if selected_raw_start < selected_raw_end:
                    selected_ref_start = overlap_ref_start + (selected_raw_start - overlap_raw_start)
                    selected_ref_end = selected_ref_start + (selected_raw_end - selected_raw_start)
                    q0 = query_position + (selected_ref_start - block_ref_start)
                    q1 = query_position + (selected_ref_end - block_ref_start)
                    _select_query(q0, q1, selected_raw_start, selected_raw_end - 1)
            query_position += length
            reference_position += length
        elif operation == 1:  # insertion relative to GRCh38
            anchor = previous_raw + 1 if previous_raw is not None else _map_ref_position_to_raw_blockwise(coord, reference_position)
            if anchor is not None and int(window_start) <= anchor < int(window_end):
                _select_query(query_position, query_position + length, int(anchor), int(anchor))
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

