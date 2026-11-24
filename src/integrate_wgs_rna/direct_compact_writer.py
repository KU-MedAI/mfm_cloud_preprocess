from __future__ import annotations

import hashlib
import csv
import json
import os
import shutil
import time
import uuid
from collections import Counter, OrderedDict
from collections.abc import Iterable
from .cpu_executor import WorkerBudget, cpu_results, elastic_cpu_enabled
from .compact.read_index import ReadIndex
from .compact.runtime import nonnegative_env_int
from .window_quality import depth_metrics, write_quality
from .compact.memory_preflight import forecast_compact
from pathlib import Path
from typing import Any

import numpy as np
import pysam
import zarr
from numcodecs import Blosc, blosc

from pre_training.preprocess_raw_readset.config import ExportConfig
from pre_training.preprocess_raw_readset.filters import rejection_reason
from pre_training.preprocess_raw_readset.shared_hap.manifest import (
    init_manifest,
    mark_manifest_complete,
    mark_manifest_failed,
    update_manifest_record,
    utc_now,
    write_json,
)
from pre_training.preprocess_raw_readset.shared_hap.records import hap_paths
from pre_training.preprocess_raw_readset.shared_hap.schema import HAP_LABELS
from pre_training.preprocess_raw_readset.shared_hap_cli import _job_schedule_key
from pre_training.preprocess_raw_readset.shared_hap_writer import (
    HAP_PROVENANCE_ATTR_KEYS,
    build_paired_hap_window_artifacts,
    prepare_hap_context,
)

from .legacy_shared_hap_converter import (
    READ_ALLOWED,
    SEQUENCE_ALLOWED,
    append_1d,
    directory_size,
    format_symbols,
    human_bytes,
    pack_4bit_acgtnp_chunk,
    write_windows_tsv,
)


# Fast-write zarr chunk preset.
#
# These values do not change the model window size or output schema. They only
# control how zarr splits arrays into compressed on-disk chunks. The current
# priority is faster preprocessing and fewer Lustre metadata operations while
# preserving the existing dataloader-facing layout.
READ_META_CHUNK = 524_288
SEGMENT_CHUNK = 1_048_576
PACKED_BASE_CHUNK = 32 << 20
WINDOW_META_CHUNK = 1024
DEPTH_CHUNK = 1 << 20
FULL_DIRECT_CONTIG_SCHEMA = "mfm_integrated_direct_shared_hap_4bit_acgtnp_depth_contig_v1"
FULL_DIRECT_PERSONAL_CONTIG_SCHEMA = "mfm_integrated_direct_personal_realigned_hap_split_4bit_acgtnp_depth_contig_v1"
DEFAULT_COMPACT_COMPRESSOR = "zstd"
DEFAULT_COMPACT_COMPRESSION_LEVEL = 5
SUPPORTED_COMPACT_COMPRESSORS = {"zstd", "lz4", "blosclz"}

PERSONAL_INPUT_PROVENANCE_SCHEMA = "stage70_personal_input_provenance_v1"
PERSONAL_FULL_VALIDATION_SCHEMA = "stage70_personal_full_validation_v1"


def _positive_env_int(name: str, default: int) -> int:
    """Validation read-ahead only: output schema/content does not change."""
    raw = os.environ.get(name, str(default))
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(f"{name} must be a positive integer: {raw}") from exc
    if value <= 0:
        raise ValueError(f"{name} must be a positive integer: {raw}")
    return value


VALIDATION_RECORD_CHUNK = _positive_env_int("MFM_VALIDATION_RECORD_CHUNK", 262_144)
VALIDATION_PACKED_BYTE_CHUNK = _positive_env_int("MFM_VALIDATION_PACKED_BYTE_CHUNK", 16 << 20)
VALIDATION_CACHE_CHUNKS = nonnegative_env_int("MFM_VALIDATION_CACHE_CHUNKS", 4)
FOURBIT_CODE_TO_ASCII = np.frombuffer(b"ACGTNP", dtype=np.uint8)
blosc.set_nthreads(max(1, int(os.environ.get("MFM_BLOSC_THREADS", "8"))))


class _ReadValidationChunkCache:
    """Bounded cache for contiguous read metadata and CIGAR Zarr slices.

    Full validation touches each read once for CIGAR integrity and then revisits
    reads while rebuilding per-window depth. Loading a complete, contiguous
    record chunk keeps those passes bounded in memory without issuing scalar
    Zarr reads for every field of every segment.
    """

    _META_NAMES = ("cigar_offset", "cigar_len", "ref_start", "ref_end", "lengths", "segment_type")

    def __init__(
        self,
        read_arrays: dict[str, Any],
        cigar_ops: Any,
        record_count: int,
        *,
        chunk_size: int = VALIDATION_RECORD_CHUNK,
        max_chunks: int = VALIDATION_CACHE_CHUNKS,
    ) -> None:
        self.read_arrays = read_arrays
        self.cigar_ops = cigar_ops
        self.record_count = int(record_count)
        self.chunk_size = max(1, int(chunk_size))
        self.max_chunks = int(max_chunks)
        if self.max_chunks < 0:
            raise ValueError("캐시 chunk 수는 음수일 수 없습니다.")
        self._cache: OrderedDict[int, dict[str, Any]] = OrderedDict()

    def _load(self, chunk_index: int) -> dict[str, Any]:
        lo = chunk_index * self.chunk_size
        hi = min(self.record_count, lo + self.chunk_size)
        metadata = {
            name: np.asarray(self.read_arrays[name][lo:hi])
            for name in self._META_NAMES
        }
        offsets = metadata["cigar_offset"].astype(np.int64, copy=False)
        lengths = metadata["cigar_len"].astype(np.int64, copy=False)
        if offsets.size:
            cigar_lo = int(offsets.min())
            cigar_hi = int(np.max(offsets + lengths))
        else:
            cigar_lo = cigar_hi = 0
        metadata.update(
            {
                "record_lo": lo,
                "cigar_lo": cigar_lo,
                "cigar_ops": np.asarray(self.cigar_ops[cigar_lo:cigar_hi], dtype=np.uint32),
            }
        )
        return metadata

    def record(self, read_id: int) -> tuple[dict[str, Any], int]:
        if read_id < 0 or read_id >= self.record_count:
            raise ValueError(f"full validation invalid read_id={read_id}")
        chunk_index = read_id // self.chunk_size
        if not self.max_chunks:
            chunk = self._load(chunk_index)
            return chunk, read_id - int(chunk["record_lo"])
        chunk = self._cache.pop(chunk_index, None)
        if chunk is None:
            chunk = self._load(chunk_index)
        self._cache[chunk_index] = chunk
        while len(self._cache) > self.max_chunks:
            self._cache.popitem(last=False)
        return chunk, read_id - int(chunk["record_lo"])

    def decode(self, read_id: int) -> tuple[list[tuple[int, int]], int, int, int, int]:
        chunk, local_index = self.record(read_id)
        cigar_offset = int(chunk["cigar_offset"][local_index])
        cigar_len = int(chunk["cigar_len"][local_index])
        if cigar_len <= 0:
            raise ValueError(f"full validation empty CIGAR for read_id={read_id}")
        op_lo = cigar_offset - int(chunk["cigar_lo"])
        ops = chunk["cigar_ops"][op_lo:op_lo + cigar_len]
        if ops.shape[0] != cigar_len:
            raise ValueError(f"full validation truncated CIGAR for read_id={read_id}")
        ref_start = int(chunk["ref_start"][local_index])
        ref_end = int(chunk["ref_end"][local_index])
        read_length = int(chunk["lengths"][local_index])
        segment_type = int(chunk["segment_type"][local_index])
        ref_pos = ref_start
        reference_consumed = 0
        query_consumed = 0
        intervals: list[tuple[int, int]] = []
        for packed_value in ops:
            packed = int(packed_value)
            op = packed & 0x0F
            length = packed >> 4
            if op not in {0, 1, 2, 3, 6, 7, 8} or length <= 0:
                raise ValueError(f"full validation invalid stored CIGAR op={op} length={length} read_id={read_id}")
            if op in {0, 1, 7, 8}:
                query_consumed += length
            if op in {0, 2, 3, 7, 8}:
                reference_consumed += length
            if op in {0, 7, 8}:
                intervals.append((ref_pos, ref_pos + length))
                ref_pos += length
            elif op in {2, 3}:
                ref_pos += length
        if query_consumed != read_length:
            raise ValueError(
                f"full validation CIGAR/query mismatch read_id={read_id}: cigar={query_consumed} stored={read_length}"
            )
        expected_reference = ref_end - ref_start
        if reference_consumed != expected_reference:
            raise ValueError(
                f"full validation CIGAR/reference mismatch read_id={read_id}: "
                f"cigar={reference_consumed} stored={expected_reference}"
            )
        return intervals, query_consumed, reference_consumed, read_length, segment_type

class DirectTiming:
    def __init__(self) -> None:
        self.steps: dict[str, float] = {}
        self.chunks: list[dict[str, Any]] = []
        self.counts: dict[str, int] = {}

    def add(self, name: str, elapsed: float) -> None:
        self.steps[name] = self.steps.get(name, 0.0) + float(elapsed)

    def add_count(self, name: str, value: int) -> None:
        self.counts[name] = self.counts.get(name, 0) + int(value)

    def step(self, name: str):
        return _DirectTimingStep(self, name)

    def add_chunk(self, chunk: dict[str, Any]) -> None:
        self.chunks.append(_json_ready(chunk))

    def to_dict(self) -> dict[str, Any]:
        return {
            "steps": {key: round(value, 6) for key, value in sorted(self.steps.items())},
            "counts": dict(sorted(self.counts.items())),
            "chunks": self.chunks,
        }


class _DirectTimingStep:
    def __init__(self, timing: DirectTiming, name: str) -> None:
        self.timing = timing
        self.name = name
        self.started = 0.0

    def __enter__(self):
        self.started = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.timing.add(self.name, time.perf_counter() - self.started)


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


def _timing_add_count(timing: dict[str, Any] | None, name: str, value: int) -> None:
    if timing is None:
        return
    counts = timing.setdefault("counts", {})
    counts[name] = int(counts.get(name, 0)) + int(value)


def _merge_timing_dict(target: DirectTiming, source: dict[str, Any]) -> None:
    # Accumulate per-chunk/per-label timing into the contig-level summary.
    # Detailed per-chunk payloads are optional, but aggregate timing is cheap
    # enough to keep for every direct compact run.
    for name, elapsed in source.get("steps", {}).items():
        target.add(str(name), float(elapsed))
    for name, value in source.get("counts", {}).items():
        target.add_count(str(name), int(value))


def _json_ready(value):
    if isinstance(value, dict):
        return {str(k): _json_ready(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json_ready(v) for v in value]
    if isinstance(value, tuple):
        return [_json_ready(v) for v in value]
    if isinstance(value, Counter):
        return {str(k): int(v) for k, v in value.items()}
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value




def _canonical_json_sha256(payload: dict[str, Any]) -> str:
    encoded = json.dumps(_json_ready(payload), sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _file_identity(path: Path, cache: dict[tuple[str, int, int], dict[str, Any]]) -> dict[str, Any]:
    """Return a content fingerprint; cache avoids rehashing shared CRAMs per contig."""
    resolved = path.resolve()
    stat = resolved.stat()
    key = (str(resolved), int(stat.st_size), int(stat.st_mtime_ns))
    cached = cache.get(key)
    if cached is not None:
        return dict(cached)
    digest = hashlib.sha256()
    with resolved.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    payload = {
        "path": str(resolved),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "sha256": digest.hexdigest(),
    }
    cache[key] = payload
    return dict(payload)


def _fasta_contig_identity(path: Path, contig: str) -> dict[str, Any]:
    """Fingerprint one FASTA contig by biological sequence, not FASTA formatting."""
    with pysam.FastaFile(str(path)) as fasta:
        if contig not in fasta.references:
            raise ValueError(f"FASTA contig missing: path={path} contig={contig} references={','.join(fasta.references)}")
        length = int(fasta.get_reference_length(contig))
        digest = hashlib.md5()
        for start in range(0, length, 1 << 20):
            digest.update(fasta.fetch(contig, start, min(length, start + (1 << 20))).upper().encode("ascii"))
    return {
        "path": str(path.resolve()),
        "contig": contig,
        "length": length,
        "sequence_md5": digest.hexdigest(),
    }


def _alignment_contig_identity(alignment: Path, reference_fasta: Path, contig: str) -> dict[str, Any]:
    """Verify the CRAM header axis against the exact personal FASTA contig."""
    with pysam.AlignmentFile(str(alignment), "rb", reference_filename=str(reference_fasta)) as handle:
        if contig not in handle.references:
            raise ValueError(f"personal alignment contig missing: path={alignment} contig={contig}")
        header = handle.header.to_dict()
        sq_rows = [row for row in header.get("SQ", []) if str(row.get("SN", "")) == contig]
        if len(sq_rows) != 1:
            raise ValueError(f"personal alignment SQ header ambiguity: path={alignment} contig={contig}")
        row = sq_rows[0]
        length = int(handle.get_reference_length(contig))
        sequence_md5 = str(row.get("M5", "")).lower()
    if not sequence_md5:
        raise ValueError(f"personal CRAM SQ header lacks M5: path={alignment} contig={contig}")
    return {
        "path": str(alignment.resolve()),
        "contig": contig,
        "length": length,
        "sequence_md5": sequence_md5,
    }


def _personal_input_provenance(job: dict[str, Any], file_cache: dict[tuple[str, int, int], dict[str, Any]]) -> dict[str, Any]:
    """Build an immutable input contract before any Stage 70 output is reused."""
    contig = str(job["contig"])
    labels: dict[str, Any] = {}
    for label in HAP_LABELS:
        hap_number = label[-1]
        alignment = Path(str(job[f"hap{hap_number}_alignment"]))
        target_fasta = _personal_alignment_reference(alignment)
        if target_fasta is None:
            raise ValueError(f"personal_realigned Stage 70 requires CRAM with target.fa: {alignment}")
        hap_fasta = Path(str(job[f"hap{hap_number}_fasta"]))
        hap_identity = _fasta_contig_identity(hap_fasta, contig)
        target_identity = _fasta_contig_identity(target_fasta, contig)
        alignment_identity = _alignment_contig_identity(alignment, target_fasta, contig)
        expected = (hap_identity["length"], hap_identity["sequence_md5"])
        observed_target = (target_identity["length"], target_identity["sequence_md5"])
        observed_header = (alignment_identity["length"], alignment_identity["sequence_md5"])
        if observed_target != expected:
            raise ValueError(
                f"personal HAP FASTA != CRAM target FASTA: {label}/{contig} hap={expected} target={observed_target}"
            )
        if observed_header != expected:
            raise ValueError(
                f"personal HAP FASTA != CRAM SQ header: {label}/{contig} hap={expected} header={observed_header}"
            )
        labels[label] = {
            "alignment": _file_identity(alignment, file_cache),
            "target_fasta": {**_file_identity(target_fasta, file_cache), **target_identity},
            "hap_fasta": {**_file_identity(hap_fasta, file_cache), **hap_identity},
            "alignment_header": alignment_identity,
        }
    window_index = _window_index_tsv(job.get("window_index_dir", ""))
    if window_index is None or not window_index.exists():
        raise FileNotFoundError(f"personal_realigned Stage 70 requires windows.tsv: {window_index}")
    payload = {
        "schema": PERSONAL_INPUT_PROVENANCE_SCHEMA,
        "validation_schema": PERSONAL_FULL_VALIDATION_SCHEMA,
        "contig": contig,
        "read_coordinate_mode": str(job["read_coordinate_mode"]),
        "validation_level": str(job["validation_level"]),
        "compact_compressor": str(job["compact_compressor"]),
        "compact_compression_level": int(job["compact_compression_level"]),
        "max_windows": None if job.get("max_windows") is None else int(job["max_windows"]),
        "export_config": dict(job["config"]),
        "window_index": _file_identity(window_index, file_cache),
        "haplotypes": labels,
    }
    payload["fingerprint_sha256"] = _canonical_json_sha256(payload)
    return payload


def _assert_resume_provenance(path: Path, expected: dict[str, Any]) -> None:
    metadata = json.loads((path / "contig_metadata.json").read_text(encoding="utf-8"))
    observed = metadata.get("input_provenance")
    if not isinstance(observed, dict):
        raise ValueError(
            f"existing compact has no input provenance and cannot be safely resumed: {path}; "
            "choose a new output path or explicitly set COMPACT_OVERWRITE=1"
        )
    if observed.get("fingerprint_sha256") != expected.get("fingerprint_sha256"):
        raise ValueError(
            f"existing compact input fingerprint differs: {path}; source CRAM/FASTA/window/config changed. "
            "choose a new output path or explicitly set COMPACT_OVERWRITE=1"
        )

def make_compact_compressor(name: str = DEFAULT_COMPACT_COMPRESSOR, level: int = DEFAULT_COMPACT_COMPRESSION_LEVEL) -> Blosc:
    """Return the Blosc compressor used for model-input compact zarr arrays.

    The output schema stores compact 4bit arrays regardless of the compressor.
    This option only changes the zarr chunk compression tradeoff:

    - zstd level 5: smaller output, slower write.
    - zstd level 1/3 or lz4: faster write, usually larger output.
    """
    normalized = (name or DEFAULT_COMPACT_COMPRESSOR).lower()
    if normalized not in SUPPORTED_COMPACT_COMPRESSORS:
        raise ValueError(f"unsupported compact compressor: {name!r}; expected one of {sorted(SUPPORTED_COMPACT_COMPRESSORS)}")
    clevel = int(level)
    if clevel < 0 or clevel > 9:
        raise ValueError(f"compact compression level must be 0..9: {level}")
    return Blosc(cname=normalized, clevel=clevel, shuffle=Blosc.BITSHUFFLE)


def _cigar_ref_match_intervals(record: dict[str, Any]) -> list[tuple[int, int]]:
    """Return absolute reference/HAP intervals covered by read bases.

    CIGAR semantics for depth:
    - M/= /X consume query and reference, so these bases contribute depth.
    - D/N consume reference only, so the read has no base at those positions.
    - I/S/H/P consume no reference depth position here. S/H are normally removed
      before storage, but skipping them keeps the helper robust.
    """
    ref_pos = int(record.get("reference_start", -1))
    ref_end = int(record.get("reference_end", ref_pos))
    if ref_pos < 0 or ref_end < ref_pos:
        return []
    ops = record.get("cigar_ops_packed", [])
    if ops is None or len(ops) == 0:
        return [(ref_pos, ref_end)] if ref_end > ref_pos else []
    intervals: list[tuple[int, int]] = []
    for packed_value in ops:
        packed = int(packed_value)
        op = packed & 0xF
        length = packed >> 4
        if length <= 0:
            continue
        if op in (0, 7, 8):  # M, =, X
            intervals.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op in (2, 3):  # D, N
            ref_pos += length
        elif op in (1, 4, 5, 6):  # I, S, H, P
            continue
        else:
            # Unknown/obsolete CIGAR op: do not invent reference coverage.
            continue
    return intervals


def _depth_arrays_for_windows(
    windows: list[Any],
    records_by_window: list[list[dict[str, Any]]],
    *,
    method: str,
    segment_types: tuple[int, ...] = (0, 1),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    """Build flat uint32 depth and per-window offsets/lengths.

    For personal_realigned data, ``method`` is CIGAR-aware and uses absolute
    HAP coordinates. For reference_projected data, exact CIGAR-to-HAP projection
    is not representable from the stored original CIGAR alone after indels in
    the HAP chain, so this writer falls back to segment span depth. v2
    production uses personal_realigned mode for depth intended for modelling.
    """
    flat_parts: list[np.ndarray] = []
    allowed_segment_types = {int(value) for value in segment_types}
    offsets = [0]
    lengths: list[int] = []
    depth_sum = 0
    depth_max = 0
    for window, records in zip(windows, records_by_window):
        # SequenceWindow.sequence is padded with "P" to the model window size,
        # while actual_length is the number of real FASTA bases before padding.
        # Depth is a model-input track and must share the same axis as the
        # stored sequence. Padding bases are stored with depth 0, while
        # CIGAR-derived coverage is accumulated only over actual_length.
        storage_length = len(str(window.sequence))
        actual_length = int(window.actual_length)
        lengths.append(storage_length)
        if storage_length <= 0:
            flat_parts.append(np.asarray([], dtype=np.uint32))
            offsets.append(offsets[-1])
            continue
        gpu_depth=os.environ.get('MFM_COVERAGE_BACKEND','cpu')=='cupy'
        starts=[];ends=[]
        diff = None if gpu_depth else np.zeros(storage_length + 1, dtype=np.int64)
        window_start = int(window.raw_start)
        actual_window_end = window_start + actual_length
        for record in records:
            if int(record.get("segment_type", 0)) not in allowed_segment_types:
                continue
            if method == "cigar_ref_match":
                intervals = _cigar_ref_match_intervals(record)
                for ref_start, ref_end in intervals:
                    start = max(int(ref_start), window_start)
                    end = min(int(ref_end), actual_window_end)
                    if end <= start:
                        continue
                    if gpu_depth:starts.append(start-window_start);ends.append(end-window_start)
                    else:
                        diff[start - window_start] += 1
                        diff[end - window_start] -= 1
            else:
                # Reference-projected fallback: uses already projected window
                # segment interval. This is kept for backwards-compatible
                # reference_projected outputs, not as the v2 preferred depth.
                start = max(int(record.get("read_start", 0)), 0)
                end = min(int(record.get("read_end", start)), actual_length)
                if end <= start:
                    continue
                if gpu_depth:starts.append(start);ends.append(end)
                else:
                    diff[start] += 1
                    diff[end] -= 1
        if gpu_depth:
            from .gpu_arrays import depth_from_intervals
            depth=depth_from_intervals(starts,ends,storage_length)
        else:
            depth = np.cumsum(diff[:-1], dtype=np.int64).astype(np.uint32, copy=False)
        if depth.size:
            depth_sum += int(depth.sum(dtype=np.uint64))
            depth_max = max(depth_max, int(depth.max()))
        flat_parts.append(depth)
        offsets.append(offsets[-1] + storage_length)
    flat = np.concatenate(flat_parts) if flat_parts else np.asarray([], dtype=np.uint32)
    return flat, np.asarray(offsets, dtype=np.int64), np.asarray(lengths, dtype=np.int32), depth_sum, depth_max


def run_direct_compact_jobs(
    *,
    records: list[dict[str, str]],
    record_tsv: Path,
    output_root: Path,
    windows_per_chunk: int,
    jobs: int,
    schedule: str,
    max_windows: int | None,
    validation_level: str,
    overwrite: bool,
    resume: bool,
    reference_fasta: Path | None,
    export_config: ExportConfig,
    sample_metadata: dict[str, Any],
    window_index_dir: str = "",
    read_coordinate_mode: str = "reference_projected",
    hap1_alignment: str = "",
    hap2_alignment: str = "",
    compact_compressor: str = DEFAULT_COMPACT_COMPRESSOR,
    compact_compression_level: int = DEFAULT_COMPACT_COMPRESSION_LEVEL,
    profile_direct_compact: bool = False,
) -> dict[str, Any]:
    run_jobs = [
        {
            "sample": row["sample"],
            "contig": row["contig"],
            "hap1_fasta": row["hap1_fasta"],
            "hap1_chain": row["hap1_chain"],
            "hap2_fasta": row["hap2_fasta"],
            "hap2_chain": row["hap2_chain"],
            "bam": row["bam"],
            "output_root": str(output_root),
            "windows_per_chunk": int(windows_per_chunk),
            "resume": bool(resume),
            "overwrite": bool(overwrite),
            "max_windows": max_windows,
            "validation_level": validation_level,
            "compact_compressor": compact_compressor,
            "compact_compression_level": int(compact_compression_level),
            "profile_direct_compact": bool(profile_direct_compact),
            "reference_fasta": str(reference_fasta) if reference_fasta else "",
            "window_index_dir": str(row.get("window_index_dir") or window_index_dir or ""),
            "read_coordinate_mode": str(row.get("read_coordinate_mode") or read_coordinate_mode or "reference_projected"),
            "hap1_alignment": str(row.get("hap1_alignment") or hap1_alignment or ""),
            "hap2_alignment": str(row.get("hap2_alignment") or hap2_alignment or ""),
            "config": export_config.to_dict(),
            "sample_metadata": dict(sample_metadata),
        }
        for row in records
    ]
    # provenance hashing은 부모 한 프로세스에서 수행한다. 이 구간도 여분
    # 예약을 반환하고 실제 contig 계산을 시작할 때 다시 자원을 요청한다.
    if elastic_cpu_enabled("MFM_COMPACT_JOBS_REQUESTED",thread_cost_accounted=True):
        WorkerBudget(1, 1, True).capacity(1, grow=False)
    memory_forecast = forecast_compact(
        run_jobs, windows_per_chunk, export_config.window_size,
        bounded_index=os.environ.get("MFM_COMPACT_INDEX_MODE", "bounded") == "bounded",
    )
    input_fingerprint_cache: dict[tuple[str, int, int], dict[str, Any]] = {}
    for job in run_jobs:
        if str(job["read_coordinate_mode"]) == "personal_realigned":
            job["input_provenance"] = _personal_input_provenance(job, input_fingerprint_cache)
    # 전체 파일 해싱 중에는 계산 CPU를 잡지 않고, 계산 직전에 원자 재확보한다.
    from .compact.resources import reserve_shared
    per_hap_ceiling=memory_forecast['worker_ceiling_bytes']
    execution_resources=reserve_shared(jobs,len(run_jobs),memory_forecast,
        hap_split=all(j['read_coordinate_mode']=='personal_realigned' for j in run_jobs))
    jobs=execution_resources['workers']
    for job in run_jobs:
        job.update(alignment_threads=execution_resources['alignment_threads'],
                   blosc_threads=execution_resources['blosc_threads'],
                   hap_jobs=execution_resources['hap_jobs'],
                   hap_cpu_cost=execution_resources['hap_cpu_cost'],
                   cpu_per_worker=execution_resources['cpu_per_worker'],
                   hap_memory_ceiling_bytes=memory_forecast.get('per_hap_address_space_bytes',per_hap_ceiling))
    forecast_path = output_root.parent.parent / "metadata/memory_preflight_compact.json"
    forecast_path.parent.mkdir(parents=True, exist_ok=True)
    write_json(forecast_path, memory_forecast)
    print(json.dumps({"memory_preflight":str(forecast_path),
                      "worker_ceiling_bytes":memory_forecast['worker_ceiling_bytes'],
                      "worker_reservation_bytes":memory_forecast['worker_reservation_bytes'],
                      "memory_limit_bytes":memory_forecast['memory_limit_bytes']}),flush=True)

    if schedule == "longest-first":
        run_jobs = sorted(run_jobs, key=_job_schedule_key)

    manifest = init_manifest(
        record_tsv=record_tsv,
        output_root=output_root,
        windows_per_chunk=windows_per_chunk,
        jobs=jobs,
        schedule=schedule,
        max_windows=max_windows,
        validation_level=validation_level,
        config=export_config.to_dict(),
    )
    manifest["output_layout"] = "direct_compact_full_shared_hap_4bit_acgtnp_depth"
    manifest["read_coordinate_mode"] = str(read_coordinate_mode or "reference_projected")
    manifest["compact_mode"] = "direct"
    manifest["compact_compressor"] = compact_compressor
    manifest["compact_compression_level"] = int(compact_compression_level)
    manifest["profile_direct_compact"] = bool(profile_direct_compact)
    manifest["execution_resources"]=execution_resources
    manifest_path = output_root / "manifest.json"
    write_json(manifest_path, manifest)

    try:
        for _, result, error in cpu_results(
            _run_direct_compact_contig, [(job,) for job in run_jobs], jobs,
            elastic=elastic_cpu_enabled("MFM_COMPACT_JOBS_REQUESTED",thread_cost_accounted=True),
            memory_per_worker_bytes=memory_forecast["worker_reservation_bytes"],
            address_space_per_worker_bytes=memory_forecast["worker_address_space_bytes"],
            base_memory_bytes=memory_forecast.get("stage_base_memory_bytes"),
            cpus_per_worker=execution_resources['cpu_per_worker'],
            mp_context=__import__('multiprocessing').get_context('spawn') if execution_resources['gpu'] else None,
            max_tasks_per_child=execution_resources['recycle_tasks'],
        ):
            if error is not None:
                raise error
            update_manifest_record(manifest, result)
            write_json(manifest_path, manifest)
            print(json.dumps(result, ensure_ascii=False), flush=True)
    except BaseException as exc:
        mark_manifest_failed(manifest, exc)
        write_json(manifest_path, manifest)
        raise

    mark_manifest_complete(manifest)
    write_json(manifest_path, manifest)
    return manifest


def summarize_direct_compact_outputs(
    *,
    sample_id: str,
    source_sample_id: str,
    dataset_role: str,
    read_type: str,
    project_id: str,
    output_sample_root: Path,
    compact_root: Path,
    contigs: list[str],
) -> dict[str, Any]:
    contig_summaries = []
    for contig in contigs:
        metadata_path = compact_root / contig / "contig_metadata.json"
        if not metadata_path.exists():
            raise FileNotFoundError(f"missing direct compact metadata: {metadata_path}")
        payload = json.loads(metadata_path.read_text(encoding="utf-8"))
        read_summary = payload.get("read_summary", {})
        hap_summaries = payload.get("hap_summaries", {})
        read_count = int(read_summary.get("read_count", read_summary.get("total_read_count", 0)))
        read_bases = int(read_summary.get("total_bases", 0))
        if read_bases == 0 and isinstance(read_summary.get("by_haplotype"), dict):
            read_bases = sum(
                int(item.get("total_bases", 0))
                for item in read_summary["by_haplotype"].values()
                if isinstance(item, dict)
            )
        target = compact_root / contig
        size_bytes = directory_size(target)
        contig_summaries.append(
            {
                "contig": contig,
                "status": payload.get("status", "direct_compact"),
                "source": payload.get("source_alignment", ""),
                "target": str(target),
                "window_rows": int(payload.get("window_rows", 0)),
                "read_count": read_count,
                "read_bases": read_bases,
                "has_depth": bool((target / "depth.zarr").exists()),
                "depth_summary": payload.get("depth_summary", {}),
                "hap_invalid_base_count": sum(
                    int(item.get("invalid_base_count", 0))
                    for item in hap_summaries.values()
                    if isinstance(item, dict)
                ),
                "hap_invalid_base_symbols": merge_count_dicts(
                    [
                        item.get("invalid_base_symbols", {})
                        for item in hap_summaries.values()
                        if isinstance(item, dict)
                    ]
                ),
                "size_bytes": size_bytes,
                "size_human": human_bytes(size_bytes),
            }
        )
    return {
        "sample_id": sample_id,
        "source_sample_id": source_sample_id or sample_id,
        "dataset_role": dataset_role,
        "read_type": read_type,
        "project_id": project_id,
        "output_sample_root": str(output_sample_root),
        "compact_root": str(compact_root),
        "encoding": "4bit_acgtnp",
        "compact_mode": "direct",
        "sequence_allowed_symbols": format_symbols(SEQUENCE_ALLOWED),
        "read_allowed_symbols": format_symbols(READ_ALLOWED),
        "hap_invalid_policy": "to_n",
        "hap_invalid_replacement": "N",
        "contigs": contig_summaries,
    }


def _run_direct_compact_contig(job: dict[str, Any]) -> dict[str, Any]:
    os.environ['MFM_COMPACT_ALIGNMENT_THREADS']=str(job.get('alignment_threads',1))
    os.environ['MFM_SAMPLE_CPU_BUDGET']=str(job.get('cpu_per_worker',1))
    blosc.set_nthreads(int(job.get('blosc_threads',os.environ.get('MFM_BLOSC_THREADS','1'))))
    started = time.time()
    timing = DirectTiming()
    config = ExportConfig(**job["config"])
    sample = str(job["sample"])
    contig = str(job["contig"])
    output_root = Path(job["output_root"])
    if str(job.get("read_coordinate_mode") or "reference_projected") == "personal_realigned":
        return _run_personal_realigned_compact_contig(job)
    profile_direct = bool(job.get("profile_direct_compact"))
    out_path = output_root / contig
    partial = output_root / f".{contig}.partial.{os.getpid()}"
    with timing.step("direct_prepare_output"):
        if partial.exists():
            shutil.rmtree(partial)
        if out_path.exists():
            if job["overwrite"]:
                shutil.rmtree(out_path)
            elif job["resume"] and _existing_direct_complete(out_path):
                metadata = json.loads((out_path / "contig_metadata.json").read_text(encoding="utf-8"))
                read_summary = metadata.get("read_summary", {})
                return {
                    "status": "skip_existing",
                    "sample": sample,
                    "contig": contig,
                    "zarr": str(out_path),
                    "shared_read_count": int(read_summary.get("read_count", 0)),
                    "elapsed_seconds": 0.0,
                }
            else:
                raise FileExistsError(f"existing direct compact output: {out_path}")

    try:
        summaries: dict[str, dict[str, Any]] = {}
        hap_window_counts: dict[str, int] = {}
        hap_contexts: dict[str, dict[str, Any]] = {}
        for label in HAP_LABELS:
            with timing.step(f"prepare_hap_context_{label}"):
                fasta, chain = hap_paths(job, label)
                context = prepare_hap_context(
                    sample=sample,
                    label=label,
                    contig=contig,
                    individual_fasta=fasta,
                    chain_path=chain,
                    config=config,
                    window_index_path=_window_index_tsv(job.get("window_index_dir", "")),
                )
                total = int(context["summary"]["coordinate_eligible_window_count_total"])
                if job.get("max_windows") is not None:
                    total = min(total, int(job["max_windows"]))
                if total <= 0:
                    raise ValueError(f"no eligible windows for {sample}/{label}/{contig}")
                summaries[label] = context["summary"]
                hap_window_counts[label] = total
                hap_contexts[label] = context

        with timing.step("init_direct_writer"):
            writer = DirectCompactChrWriter(
                partial,
                final_target=out_path,
                sample=sample,
                contig=contig,
                hap_window_counts=hap_window_counts,
                window_size=config.window_size,
                source_alignment=Path(job["bam"]),
                reference_fasta=Path(job["reference_fasta"]) if job.get("reference_fasta") else None,
                sample_metadata=job.get("sample_metadata", {}),
                compact_compressor=str(job.get("compact_compressor") or DEFAULT_COMPACT_COMPRESSOR),
                compact_compression_level=int(job.get("compact_compression_level") or DEFAULT_COMPACT_COMPRESSION_LEVEL),
                attrs={
                    "source_bam_path": str(job["bam"]),
                    "reference_fasta": str(job["reference_fasta"]) if job.get("reference_fasta") else "",
                    "window_index_dir": str(job.get("window_index_dir", "")),
                    "source_hap1_fasta": str(job["hap1_fasta"]),
                    "source_hap1_chain": str(job["hap1_chain"]),
                    "read_coordinate_mode": str(job.get("read_coordinate_mode") or "reference_projected"),
                    "source_hap2_fasta": str(job["hap2_fasta"]),
                    "source_hap2_chain": str(job["hap2_chain"]),
                    "config": config.to_dict(),
                    "window_summary": summaries,
                    "builder": "direct_compact_paired_hap_union_bam_fetch_hp_ps",
                },
            )
        chunk_size = int(job["windows_per_chunk"])
        max_total = max(hap_window_counts.values())
        with timing.step("process_chunks_total"):
            for chunk_index, start_index in enumerate(range(0, max_total, chunk_size)):
                chunk_timing: dict[str, Any] | None = None
                if profile_direct:
                    chunk_timing = {
                        "chunk_index": chunk_index,
                        "window_start_index": int(start_index),
                        "requested_windows": int(min(chunk_size, max_total - start_index)),
                        "append_by_label": [],
                    }
                chunk_started = time.perf_counter()
                artifact_timing: dict[str, Any] = chunk_timing if chunk_timing is not None else {}
                with timing.step("chunk_build_artifacts_total"):
                    artifacts_by_label = build_paired_hap_window_artifacts(
                        sample=sample,
                        contig=contig,
                        bam_path=Path(job["bam"]),
                        hap_contexts=hap_contexts,
                        reference_fasta=Path(job["reference_fasta"]) if job.get("reference_fasta") else None,
                        max_windows=min(chunk_size, max_total - start_index),
                        window_start_index=start_index,
                        config=config,
                        timing=artifact_timing,
                    )
                _merge_timing_dict(timing, artifact_timing)
                if chunk_timing is not None:
                    chunk_timing["window_count_by_hap"] = {
                        label: int(artifacts.get("window_count", 0))
                        for label, artifacts in artifacts_by_label.items()
                    }
                    chunk_timing["segment_count_by_hap"] = {
                        label: int(artifacts.get("segment_count", 0))
                        for label, artifacts in artifacts_by_label.items()
                    }
                for label, artifacts in artifacts_by_label.items():
                    if writer.hap_window_cursor[label] < hap_window_counts[label]:
                        label_timing: dict[str, Any] = {"label": label}
                        with timing.step("chunk_append_hap_total"):
                            writer.append_hap_chunk(label, artifacts, timing=label_timing)
                        _merge_timing_dict(timing, label_timing)
                        if chunk_timing is not None:
                            chunk_timing["append_by_label"].append(label_timing)
                if chunk_timing is not None:
                    chunk_timing["elapsed_seconds"] = round(time.perf_counter() - chunk_started, 6)
                    timing.add_chunk(chunk_timing)
        with timing.step("finalize_contig"):
            contig_summary = writer.finalize(validation_level=job["validation_level"], timing=timing)
        with timing.step("rename_partial"):
            partial.rename(out_path)
        write_json(out_path / "contig_timing.json", timing.to_dict())
        elapsed = round(time.time() - started, 3)
        result = {
            "status": "pass",
            "sample": sample,
            "contig": contig,
            "zarr": str(out_path),
            "shared_read_count": contig_summary.get("read_summary", {}).get("read_count", 0),
            "haplotypes": {
                label: {
                    "windows": int(item.get("window_count", 0)),
                    "segments": int(item.get("segment_count", 0)),
                }
                for label, item in contig_summary.get("hap_summaries", {}).items()
            },
            "read_haplotype_counts": contig_summary.get("read_summary", {}).get("hp_counts", {}),
            "read_phase_set": contig_summary.get("read_summary", {}).get("ps_summary", {}),
            "elapsed_seconds": elapsed,
            "direct_timing_summary": timing.to_dict().get("steps", {}),
            "direct_timing_counts": timing.to_dict().get("counts", {}),
        }
        return result
    except BaseException:
        if partial.exists():
            shutil.rmtree(partial)
        raise


def _window_index_tsv(value: str | Path | None) -> Path | None:
    if value is None or not str(value):
        return None
    path = Path(value)
    if path.is_dir():
        return path / "windows.tsv"
    return path

def _existing_direct_complete(path: Path) -> bool:
    return (
        path.exists()
        and (path / "contig_metadata.json").exists()
        and (path / "hap_sequence.zarr").exists()
        and (path / "readset.zarr").exists()
        and (path / "depth.zarr").exists()
        and (path / "windows.tsv").exists()
    )


def _existing_personal_realigned_complete(path: Path) -> bool:
    return (
        path.exists()
        and (path / "contig_metadata.json").exists()
        and (path / "hap_sequence.zarr").exists()
        and (path / "readset.HAP1.zarr").exists()
        and (path / "readset.HAP2.zarr").exists()
        and (path / "depth.zarr").exists()
        and (path / "windows.tsv").exists()
    )


def _personal_alignment_reference(alignment: Path) -> Path | None:
    if alignment.suffix.lower() != ".cram":
        return None
    ref = alignment.parent / "target.fa"
    if not ref.exists():
        raise FileNotFoundError(f"personal realigned CRAM requires target.fa next to CRAM: {ref}")
    return ref


def _run_personal_realigned_compact_contig(job: dict[str, Any]) -> dict[str, Any]:
    started = time.time()
    timing = DirectTiming()
    config = ExportConfig(**job["config"])
    sample = str(job["sample"])
    contig = str(job["contig"])
    output_root = Path(job["output_root"])
    out_path = output_root / contig
    partial = output_root / f".{contig}.personal.partial.{os.getpid()}.{uuid.uuid4().hex[:8]}"
    profile_direct = bool(job.get("profile_direct_compact"))

    input_provenance = job.get("input_provenance")
    if not isinstance(input_provenance, dict) or not input_provenance.get("fingerprint_sha256"):
        raise ValueError(f"personal_realigned Stage 70 missing input provenance: sample={sample} contig={contig}")
    with timing.step("direct_prepare_output"):
        if partial.exists():
            shutil.rmtree(partial)
        if out_path.exists():
            if job["overwrite"]:
                shutil.rmtree(out_path)
            elif job["resume"] and _existing_personal_realigned_complete(out_path):
                _assert_resume_provenance(out_path, input_provenance)
                metadata = json.loads((out_path / "contig_metadata.json").read_text(encoding="utf-8"))
                read_summary = metadata.get("read_summary", {})
                return {
                    "status": "skip_existing",
                    "sample": sample,
                    "contig": contig,
                    "zarr": str(out_path),
                    "readset_layout": "hap_split",
                    "shared_read_count": int(read_summary.get("total_read_count", 0)),
                    "elapsed_seconds": 0.0,
                }
            else:
                raise FileExistsError(f"existing personal realigned compact output: {out_path}")
        partial.mkdir(parents=True, exist_ok=True)

    try:
        alignments: dict[str, Path] = {}
        alignment_refs: dict[str, Path | None] = {}
        for label in HAP_LABELS:
            hap_number = label[-1]
            alignment_value = str(job.get(f"hap{hap_number}_alignment") or "")
            if not alignment_value:
                raise ValueError(f"personal_realigned mode requires hap{hap_number}_alignment for {sample}/{contig}")
            alignment_path = Path(alignment_value)
            if not alignment_path.exists():
                raise FileNotFoundError(f"missing {label} personal alignment: {alignment_path}")
            alignments[label] = alignment_path
            alignment_refs[label] = _personal_alignment_reference(alignment_path)
        per_hap_metadata: dict[str, dict[str, Any]] = {}
        from .compact.personal_hap import write_hap,parallel_haps
        tasks=[(job,label,str(partial/f'.{label}.single')) for label in HAP_LABELS]
        if int(job.get('hap_jobs',1))>1:
            results=parallel_haps(tasks)
        else:
            results=((i,write_hap(*task),None) for i,task in enumerate(tasks))
        for _,result,error in results:
            if error is not None:raise error
            per_hap_metadata[result['label']]=result['metadata']
            _merge_timing_dict(timing,result['timing'])
            for chunk in result['timing'].get('chunks',[]):timing.add_chunk(chunk)
        for label in HAP_LABELS:
            single_root=partial/f'.{label}.single'
            if not (partial / "hap_sequence.zarr").exists():
                shutil.move(str(single_root / "hap_sequence.zarr"), str(partial / "hap_sequence.zarr"))
            else:
                shutil.move(str(single_root / "hap_sequence.zarr" / label), str(partial / "hap_sequence.zarr" / label))
            if not (partial / "depth.zarr").exists():
                shutil.move(str(single_root / "depth.zarr"), str(partial / "depth.zarr"))
            else:
                shutil.move(str(single_root / "depth.zarr" / label), str(partial / "depth.zarr" / label))
            if (single_root / "supplementary_depth.zarr").exists():
                if not (partial / "supplementary_depth.zarr").exists():
                    shutil.move(str(single_root / "supplementary_depth.zarr"), str(partial / "supplementary_depth.zarr"))
                else:
                    shutil.move(str(single_root / "supplementary_depth.zarr" / label), str(partial / "supplementary_depth.zarr" / label))

            shutil.move(str(single_root / "readset.zarr"), str(partial / f"readset.{label}.zarr"))
            quality = single_root / "window_quality.tsv"
            if quality.exists():
                with quality.open(newline="") as handle:
                    per_hap_metadata[label]["window_quality_rows"] = list(csv.DictReader(handle, delimiter="\t"))
            shutil.rmtree(single_root)

        hap_seq_root = zarr.open_group(str(partial / "hap_sequence.zarr"), mode="a")
        hap_seq_root.attrs.update(
            {
                "schema": "mfm_hap_sequence_4bit_acgtnp_direct_v0",
                "read_coordinate_mode": "personal_realigned",
                "readset_layout": "hap_split",
                "hap_labels": list(HAP_LABELS),
                "sample": sample,
                "contig": contig,
                "window_size": config.window_size,
                "updated_at": utc_now(),
            }
        )

        depth_root = zarr.open_group(str(partial / "depth.zarr"), mode="a")
        depth_root.attrs.update(
            {
                "schema": "mfm_depth_uint32_direct_v1",
                "read_coordinate_mode": "personal_realigned",
                "readset_layout": "hap_split",
                "depth_method": "cigar_ref_match",
                "sample": sample,
                "contig": contig,
                "window_size": config.window_size,
                "hap_labels": list(HAP_LABELS),
                "updated_at": utc_now(),
            }
        )
        depth_root.attrs["depth_role"] = "model_primary"
        depth_root.attrs["input_segment_types"] = list(config.depth_segment_types)
        if (partial / "supplementary_depth.zarr").exists():
            supplementary_root = zarr.open_group(str(partial / "supplementary_depth.zarr"), mode="a")
            supplementary_root.attrs.update(
                {
                    "schema": "mfm_supplementary_depth_uint32_direct_v1",
                    "read_coordinate_mode": "personal_realigned",
                    "readset_layout": "hap_split",
                    "depth_role": "supplementary_evidence_only",
                    "input_segment_types": list(config.supplementary_depth_segment_types),
                    "hap_labels": list(HAP_LABELS),
                }
            )


        hap_summaries = {label: per_hap_metadata[label]["hap_summaries"][label] for label in HAP_LABELS}
        read_summaries = {label: per_hap_metadata[label]["read_summary"] for label in HAP_LABELS}
        depth_summaries = {
            label: per_hap_metadata[label].get("depth_summary", {}).get("by_haplotype", {}).get(label, {})
            for label in HAP_LABELS
        }
        supplementary_depth_summaries = (
            {
                label: per_hap_metadata[label].get("depth_summary", {}).get("supplementary_evidence", {}).get("by_haplotype", {}).get(label, {})
                for label in HAP_LABELS
            }
            if config.supplementary_depth_segment_types
            else None
        )
        filter_qc_payload = {
            "schema": "stage70_alignment_filter_qc_v1",
            "read_coordinate_mode": "personal_realigned",
            "filter_preset": config.filter_preset,
            "by_haplotype": {label: per_hap_metadata[label].get("filter_qc", {}) for label in HAP_LABELS},
        }
        windows = []
        for label, summary in hap_summaries.items():
            windows.extend(make_window_rows(contig, label, summary))
        write_windows_tsv(partial / "windows.tsv", windows)
        write_quality(partial / "window_quality.tsv", [row for label in HAP_LABELS
                      for row in per_hap_metadata[label].pop("window_quality_rows", [])])
        total_read_count = sum(int(item.get("read_count", 0)) for item in read_summaries.values())
        payload = {
            "schema": FULL_DIRECT_PERSONAL_CONTIG_SCHEMA,
            "status": "direct_compact",
            "read_coordinate_mode": "personal_realigned",
            "readset_layout": "hap_split",
            "sample_id": sample,
            "source_sample_id": job.get("sample_metadata", {}).get("source_sample_id", sample),
            "dataset_role": job.get("sample_metadata", {}).get("dataset_role", ""),
            "read_type": job.get("sample_metadata", {}).get("read_type", ""),
            "project_id": job.get("sample_metadata", {}).get("project_id", ""),
            "layout_name": job.get("sample_metadata", {}).get("layout_name", ""),
            "contig": contig,
            "source_alignment": {label: str(alignments[label]) for label in HAP_LABELS},
            "reference_fasta": {label: str(alignment_refs[label]) if alignment_refs[label] else "" for label in HAP_LABELS},
            "target": str(out_path),
            "window_rows": len(windows),
            "hap_summaries": hap_summaries,
            "read_summary": {"layout": "hap_split", "total_read_count": total_read_count, "by_haplotype": read_summaries},
            "depth_summary": {
                "layout": "hap_split",
                "dtype": "uint32",
                "method": "cigar_ref_match",
                "by_haplotype": depth_summaries,
                "input_segment_types": list(config.depth_segment_types),
                "supplementary_evidence": None if supplementary_depth_summaries is None else {"input_segment_types": list(config.supplementary_depth_segment_types), "by_haplotype": supplementary_depth_summaries},
            },
            "filter_policy": {
                "preset": config.filter_preset,
                "depth_segment_types": list(config.depth_segment_types),
                "supplementary_depth_segment_types": list(config.supplementary_depth_segment_types),
                "model_segment_types": list(config.model_segment_types),
            },
            "filter_qc": filter_qc_payload,
            "encoding": "4bit_acgtnp",
            "hap_invalid_policy": "to_n",
            "hap_invalid_replacement": "N",
            "input_provenance": input_provenance,
            "validation": {
                "schema": PERSONAL_FULL_VALIDATION_SCHEMA,
                "level": job["validation_level"],
                "by_haplotype": {label: per_hap_metadata[label].get("validation", {}) for label in HAP_LABELS},
            },
            "validation_level": job["validation_level"],
            "completed_at": utc_now(),
            "direct_timing": timing.to_dict(),
        }
        write_json(partial / "filter_qc.json", filter_qc_payload)
        write_json(partial / "contig_metadata.json", payload)
        write_json(partial / "contig_timing.json", timing.to_dict())
        partial.rename(out_path)
        elapsed = round(time.time() - started, 3)
        return {
            "status": "pass",
            "sample": sample,
            "contig": contig,
            "zarr": str(out_path),
            "readset_layout": "hap_split",
            "shared_read_count": total_read_count,
            "haplotypes": {label: {"windows": int(hap_summaries[label].get("window_count", 0)), "segments": int(hap_summaries[label].get("segment_count", 0)), "reads": int(read_summaries[label].get("read_count", 0))} for label in HAP_LABELS},
            "read_haplotype_counts": {label: read_summaries[label].get("hp_counts", {}) for label in HAP_LABELS},
            "read_phase_set": {label: read_summaries[label].get("ps_summary", {}) for label in HAP_LABELS},
            "depth_summary": depth_summaries,
            "elapsed_seconds": elapsed,
            "direct_timing_summary": timing.steps,
            "direct_timing_counts": timing.counts,
        }
    except Exception as exc:
        if partial.exists():
            write_json(partial / 'failure.json', {'status': 'failed', 'error': f'{type(exc).__name__}: {exc}'})
        raise


def alignment_filter_qc(
    *, alignment: Path, reference_fasta: Path | None, contig: str, config_values: dict[str, Any]
) -> dict[str, Any]:
    """Audit the exact record-level Stage 70 gate once per source contig.

    The audit deliberately retains source `@PG`: primary/supplementary/duplicate
    flags are provenance-dependent, so downstream users can see which program
    produced the supplied alignment instead of inferring that policy from depth.
    """
    config = ExportConfig(**config_values)
    config.validate()
    payload: dict[str, Any] = {
        "schema": "stage70_alignment_filter_qc_v1",
        "alignment": str(alignment),
        "reference_fasta": str(reference_fasta) if reference_fasta else "",
        "requested_contig": contig,
        "filter_preset": config.filter_preset,
        "resolved_filter": {
            "min_mapping_quality": config.min_mapping_quality,
            "exclude_secondary": config.exclude_secondary,
            "exclude_supplementary": config.exclude_supplementary,
            "exclude_duplicate": config.exclude_duplicate,
            "exclude_qcfail": config.exclude_qcfail,
            "require_proper_pair": config.require_proper_pair,
        },
    }
    if not alignment.exists():
        return {**payload, "status": "not_available", "reason": "source_alignment_missing"}
    counts: Counter[str] = Counter()
    try:
        with pysam.AlignmentFile(
            str(alignment), "rb", reference_filename=str(reference_fasta) if reference_fasta else None
        ) as handle:
            resolved_contig = contig if contig in handle.references else contig.removeprefix("chr")
            if resolved_contig not in handle.references:
                return {**payload, "status": "not_available", "reason": "contig_missing"}
            header = handle.header.to_dict()
            for read in handle.fetch(resolved_contig):
                counts["total_records"] += 1
                for flag_name, value in (
                    ("unmapped", read.is_unmapped),
                    ("secondary", read.is_secondary),
                    ("supplementary", read.is_supplementary),
                    ("duplicate", read.is_duplicate),
                    ("qcfail", read.is_qcfail),
                ):
                    if value:
                        counts[f"flag_{flag_name}"] += 1
                if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                    counts["flag_primary"] += 1
                if int(read.mapping_quality) < int(config.min_mapping_quality):
                    counts["flag_low_mapq"] += 1
                reason = rejection_reason(read, config)
                if reason is None:
                    counts["retained_records"] += 1
                else:
                    counts[f"excluded_{reason}"] += 1
    except (OSError, ValueError) as exc:
        return {**payload, "status": "not_available", "reason": type(exc).__name__}
    return {
        **payload,
        "status": "pass",
        "resolved_contig": resolved_contig,
        "header_program_groups": header.get("PG", []),
        "counts": {key: int(value) for key, value in sorted(counts.items())},
    }


class DirectCompactChrWriter:
    def __init__(
        self,
        path: Path,
        *,
        final_target: Path,
        sample: str,
        contig: str,
        hap_window_counts: dict[str, int],
        window_size: int,
        source_alignment: Path,
        reference_fasta: Path | None,
        sample_metadata: dict[str, Any],
        compact_compressor: str = DEFAULT_COMPACT_COMPRESSOR,
        compact_compression_level: int = DEFAULT_COMPACT_COMPRESSION_LEVEL,
        attrs: dict[str, Any] | None = None,
    ) -> None:
        self.path = Path(path)
        self.final_target = Path(final_target)
        self.sample = sample
        self.contig = contig
        self.window_size = int(window_size)
        self.source_alignment = source_alignment
        self.reference_fasta = reference_fasta
        self.sample_metadata = dict(sample_metadata)
        self.compact_compressor = compact_compressor
        self.compact_compression_level = int(compact_compression_level)
        self.hap_window_counts = {label: int(count) for label, count in hap_window_counts.items()}
        self.hap_window_cursor = {label: 0 for label in self.hap_window_counts}
        self.hap_segment_cursor = {label: 0 for label in self.hap_window_counts}
        self.hap_window_content_hash: dict[str, Any] = {label: hashlib.sha256() for label in self.hap_window_counts}
        self.hap_sequence_offsets: dict[str, list[int]] = {label: [0] for label in self.hap_window_counts}
        self.hap_raw_starts: dict[str, list[int]] = {label: [] for label in self.hap_window_counts}
        self.hap_actual_lengths: dict[str, list[int]] = {label: [] for label in self.hap_window_counts}
        self.hap_total_bases: dict[str, int] = {label: 0 for label in self.hap_window_counts}
        self.hap_segment_offsets: dict[str, list[int]] = {label: [0] for label in self.hap_window_counts}
        self.hap_counts: dict[str, Counter[str]] = {label: Counter() for label in self.hap_window_counts}
        self.hap_coverage: dict[str, list[float]] = {label: [] for label in self.hap_window_counts}
        self.read_index = ReadIndex()
        self.read_key_to_index = self.read_index.keys
        self.read_haplotype_by_index = self.read_index.haplotypes
        self.read_phase_set_by_index = self.read_index.phase_sets
        self.read_lengths_by_index = self.read_index.lengths
        self.read_base_offset = 0
        self.read_hp_counts: Counter[int] = Counter()
        self.read_ps_counts: Counter[int] = Counter()
        self.read_length_sum = 0
        self.read_length_min: int | None = None
        self.read_length_max: int | None = None
        self.haplotype_conflict_count = 0
        self.phase_set_conflict_count = 0
        self.cigar_op_count = 0
        self.hap_invalid_counts: dict[str, Counter[str]] = {label: Counter() for label in self.hap_window_counts}
        self.read_coordinate_mode = str((attrs or {}).get("read_coordinate_mode") or "reference_projected")
        self.depth_method = "cigar_ref_match" if self.read_coordinate_mode == "personal_realigned" else "segment_span_projection_fallback"
        self.depth_offsets: dict[str, list[int]] = {label: [0] for label in self.hap_window_counts}
        policy_config = dict((attrs or {}).get("config") or {})
        self.filter_preset = str(policy_config.get("filter_preset") or "legacy_all_retained_v0")
        self.depth_segment_types = tuple(int(value) for value in policy_config.get("depth_segment_types", (0, 1)))
        self.supplementary_depth_segment_types = tuple(
            int(value) for value in policy_config.get("supplementary_depth_segment_types", ())
        )
        self.model_segment_types = tuple(int(value) for value in policy_config.get("model_segment_types", (0, 1)))
        if set(self.depth_segment_types) & set(self.supplementary_depth_segment_types):
            raise ValueError("Stage 70 primary/supplementary depth segment types overlap")
        self.has_supplementary_depth = bool(self.supplementary_depth_segment_types)

        self.depth_total_bases: dict[str, int] = {label: 0 for label in self.hap_window_counts}
        self.filter_qc = alignment_filter_qc(alignment=self.source_alignment, reference_fasta=self.reference_fasta, contig=self.contig, config_values=policy_config)
        self.depth_sum: dict[str, int] = {label: 0 for label in self.hap_window_counts}
        self.depth_max: dict[str, int] = {label: 0 for label in self.hap_window_counts}
        # HAP provenance attributes are chunk-invariant. Write them once per HAP
        self.supplementary_depth_offsets: dict[str, list[int]] = {label: [0] for label in self.hap_window_counts}
        self.supplementary_depth_total_bases: dict[str, int] = {label: 0 for label in self.hap_window_counts}
        self.supplementary_depth_sum: dict[str, int] = {label: 0 for label in self.hap_window_counts}
        self.supplementary_depth_max: dict[str, int] = {label: 0 for label in self.hap_window_counts}

        # to avoid repeated zarr metadata reads/checks on every append chunk.
        self.hap_provenance_written: set[str] = set()

        comp = make_compact_compressor(self.compact_compressor, self.compact_compression_level)
        self.comp = comp
        self.target_root = self.path
        self.hap_seq_root = zarr.open_group(str(self.path / "hap_sequence.zarr"), mode="w")
        self.readset_root = zarr.open_group(str(self.path / "readset.zarr"), mode="w")
        self.depth_root = zarr.open_group(str(self.path / "depth.zarr"), mode="w")
        self.supplementary_depth_root = zarr.open_group(str(self.path / "supplementary_depth.zarr"), mode="w") if self.has_supplementary_depth else None
        self.hap_seq_root.attrs.update(
            {
                "schema": "mfm_hap_sequence_4bit_acgtnp_direct_v0",
                "encoding": "4bit_acgtnp",
                "symbols": format_symbols(SEQUENCE_ALLOWED),
                "sample": sample,
                "contig": contig,
                "window_size": self.window_size,
                "hap_labels": list(self.hap_window_counts),
                "created_at": utc_now(),
                "compressor": self.compact_compressor,
                "compression_level": self.compact_compression_level,
                **(attrs or {}),
            }
        )
        self.depth_root.attrs.update(
            {
                "schema": "mfm_depth_uint32_direct_v1",
                "dtype": "uint32",
                "sample": sample,
                "contig": contig,
                "window_size": self.window_size,
                "hap_labels": list(self.hap_window_counts),
                "read_coordinate_mode": self.read_coordinate_mode,
                "depth_method": self.depth_method,
                "created_at": utc_now(),
                "compressor": self.compact_compressor,
                "compression_level": self.compact_compression_level,
                **(attrs or {}),
            }
        )
        self.depth_root.attrs["depth_role"] = "model_primary"
        self.depth_root.attrs["input_segment_types"] = list(self.depth_segment_types)
        if self.supplementary_depth_root is not None:
            self.supplementary_depth_root.attrs.update(
                {
                    "schema": "mfm_supplementary_depth_uint32_direct_v1",
                    "dtype": "uint32",
                    "sample": sample,
                    "contig": contig,
                    "window_size": self.window_size,
                    "hap_labels": list(self.hap_window_counts),
                    "depth_role": "supplementary_evidence_only",
                    "input_segment_types": list(self.supplementary_depth_segment_types),
                }
            )
        self.readset_root.attrs.update(
            {
                "schema": "mfm_readset_4bit_acgtnp_direct_v0",
                "encoding": "4bit_acgtnp",
                "symbols": format_symbols(READ_ALLOWED),
                "sample": sample,
                "contig": contig,
                "window_size": self.window_size,
                "read_sequence_storage": "shared_clipping_removed_mapped_query_sequence",
                "read_haplotype_storage": "BAM_HP_tag_uint8_zero_for_missing_or_conflict",
                "read_phase_set_storage": "BAM_PS_tag_int64_zero_for_missing_or_conflict",
                "read_cigar_storage": "uint32_packed_htslib_cigar_ops_clipping_removed_length_shift_4_or_op",
                "read_reference_storage": "bam_reference_0based_halfopen_per_shared_alignment_segment",
                "created_at": utc_now(),
                "compressor": self.compact_compressor,
                "compression_level": self.compact_compression_level,
                **(attrs or {}),
            }
        )
        self.readset_root.attrs["model_segment_types"] = list(self.model_segment_types)
        self.readset_root.attrs["model_segment_policy"] = "primary_only" if self.model_segment_types == (0,) else "explicit_segment_types"

        reads = self.readset_root.create_group("reads")
        reads.attrs.update(
            {
                "encoding": "4bit_acgtnp",
                "symbols": format_symbols(READ_ALLOWED),
                "kind": "read_sequence",
                "invalid_base_policy": "fail",
                "bases_per_byte": 2,
            }
        )
        self.read_arrays = {
            "bases_packed": reads.empty("bases_packed", shape=(0,), dtype=np.uint8, chunks=(PACKED_BASE_CHUNK,), compressor=comp),
            "offsets": reads.empty("offsets", shape=(1,), dtype=np.int64, chunks=(READ_META_CHUNK,), compressor=comp),
            "lengths": reads.empty("lengths", shape=(0,), dtype=np.int32, chunks=(READ_META_CHUNK,), compressor=comp),
            "hp": reads.empty("hp", shape=(0,), dtype=np.uint8, chunks=(READ_META_CHUNK,), compressor=comp),
            "ps": reads.empty("ps", shape=(0,), dtype=np.int64, chunks=(READ_META_CHUNK,), compressor=comp),
            "ref_start": reads.empty("ref_start", shape=(0,), dtype=np.int64, chunks=(READ_META_CHUNK,), compressor=comp),
            "ref_end": reads.empty("ref_end", shape=(0,), dtype=np.int64, chunks=(READ_META_CHUNK,), compressor=comp),
            "segment_type": reads.empty("segment_type", shape=(0,), dtype=np.int8, chunks=(READ_META_CHUNK,), compressor=comp),
            "cigar_offset": reads.empty("cigar_offset", shape=(0,), dtype=np.int64, chunks=(READ_META_CHUNK,), compressor=comp),
            "cigar_len": reads.empty("cigar_len", shape=(0,), dtype=np.int32, chunks=(READ_META_CHUNK,), compressor=comp),
        }
        self.read_arrays["offsets"][0] = 0
        cigar = self.readset_root.create_group("cigar")
        self.cigar_ops = cigar.empty("ops_packed", shape=(0,), dtype=np.uint32, chunks=(SEGMENT_CHUNK,), compressor=comp)
        self.segment_groups: dict[str, zarr.Group] = {}
        self.segment_arrays: dict[str, dict[str, zarr.Array]] = {}
        self.hap_groups: dict[str, zarr.Group] = {}
        self.hap_arrays: dict[str, dict[str, zarr.Array]] = {}
        self.depth_groups: dict[str, zarr.Group] = {}
        self.depth_arrays: dict[str, dict[str, zarr.Array]] = {}
        self.supplementary_depth_groups: dict[str, zarr.Group] = {}
        self.supplementary_depth_arrays: dict[str, dict[str, zarr.Array]] = {}
        segments = self.readset_root.create_group("segments")
        for label, total_windows in self.hap_window_counts.items():
            hap_group = self.hap_seq_root.create_group(label)
            hap_group.attrs.update(
                {
                    "encoding": "4bit_acgtnp",
                    "symbols": format_symbols(SEQUENCE_ALLOWED),
                    "kind": "hap_window_sequence",
                    "invalid_base_policy": "to_n",
                    "invalid_base_replacement": "N",
                    "bases_per_byte": 2,
                }
            )
            self.hap_groups[label] = hap_group
            self.hap_arrays[label] = {
                "bases_packed": hap_group.empty("bases_packed", shape=(0,), dtype=np.uint8, chunks=(PACKED_BASE_CHUNK,), compressor=comp),
                "offsets": hap_group.empty("offsets", shape=(1,), dtype=np.int64, chunks=(WINDOW_META_CHUNK,), compressor=comp),
                "lengths": hap_group.empty("lengths", shape=(0,), dtype=np.int32, chunks=(WINDOW_META_CHUNK,), compressor=comp),
                "raw_start": hap_group.empty("raw_start", shape=(total_windows,), dtype=np.int64, chunks=(min(total_windows + 1, WINDOW_META_CHUNK),), compressor=comp),
                "actual_length": hap_group.empty("actual_length", shape=(total_windows,), dtype=np.int64, chunks=(min(total_windows + 1, WINDOW_META_CHUNK),), compressor=comp),
            }
            self.hap_arrays[label]["offsets"][0] = 0

            depth_group = self.depth_root.create_group(label)
            depth_group.attrs.update(
                {
                    "kind": "per_base_depth",
                    "dtype": "uint32",
                    "coordinate_system": f"{label}_individual_fasta",
                    "window_order": "same_as_hap_sequence",
                    "depth_method": self.depth_method,
                }
            )
            self.depth_groups[label] = depth_group
            self.depth_arrays[label] = {
                "depth": depth_group.empty("depth", shape=(0,), dtype=np.uint32, chunks=(DEPTH_CHUNK,), compressor=comp),
                "offsets": depth_group.empty("offsets", shape=(1,), dtype=np.int64, chunks=(WINDOW_META_CHUNK,), compressor=comp),
                "lengths": depth_group.empty("lengths", shape=(0,), dtype=np.int32, chunks=(WINDOW_META_CHUNK,), compressor=comp),
            }
            self.depth_arrays[label]["offsets"][0] = 0
            if self.supplementary_depth_root is not None:
                supplementary_group = self.supplementary_depth_root.create_group(label)
                supplementary_group.attrs["input_segment_types"] = list(self.supplementary_depth_segment_types)
                self.supplementary_depth_groups[label] = supplementary_group
                self.supplementary_depth_arrays[label] = {
                    "depth": supplementary_group.empty("depth", shape=(0,), dtype=np.uint32, chunks=(DEPTH_CHUNK,), compressor=comp),
                    "offsets": supplementary_group.empty("offsets", shape=(1,), dtype=np.int64, chunks=(WINDOW_META_CHUNK,), compressor=comp),
                    "lengths": supplementary_group.empty("lengths", shape=(0,), dtype=np.int32, chunks=(WINDOW_META_CHUNK,), compressor=comp),
                }
                self.supplementary_depth_arrays[label]["offsets"][0] = 0


            seg_group = segments.create_group(label)
            self.segment_groups[label] = seg_group
            self.segment_arrays[label] = {
                "read_id": seg_group.empty("read_id", shape=(0,), dtype=np.int32, chunks=(SEGMENT_CHUNK,), compressor=comp),
                "hap_start": seg_group.empty("hap_start", shape=(0,), dtype=np.int64, chunks=(SEGMENT_CHUNK,), compressor=comp),
                "hap_end": seg_group.empty("hap_end", shape=(0,), dtype=np.int64, chunks=(SEGMENT_CHUNK,), compressor=comp),
                "query_start": seg_group.empty("query_start", shape=(0,), dtype=np.int32, chunks=(SEGMENT_CHUNK,), compressor=comp),
                "query_end": seg_group.empty("query_end", shape=(0,), dtype=np.int32, chunks=(SEGMENT_CHUNK,), compressor=comp),
                "strand": seg_group.empty("strand", shape=(0,), dtype=np.int8, chunks=(SEGMENT_CHUNK,), compressor=comp),
                "window_read_offsets": seg_group.empty("window_read_offsets", shape=(total_windows + 1,), dtype=np.int64, chunks=(min(total_windows + 1, WINDOW_META_CHUNK),), compressor=comp),
            }
            self.segment_arrays[label]["window_read_offsets"][0] = 0

    @staticmethod
    def _decrement_counter(counter: Counter[int], value: int) -> None:
        counter[int(value)] -= 1
        if counter[int(value)] <= 0:
            del counter[int(value)]

    def _replace_read_hp(self, index: int, value: int) -> None:
        old_value = int(self.read_haplotype_by_index[index])
        new_value = int(value)
        if old_value == new_value:
            return
        self._decrement_counter(self.read_hp_counts, old_value)
        self.read_hp_counts[new_value] += 1
        self.read_haplotype_by_index[index] = new_value
        self.read_arrays["hp"][index] = new_value

    def _replace_read_ps(self, index: int, value: int) -> None:
        old_value = int(self.read_phase_set_by_index[index])
        new_value = int(value)
        if old_value == new_value:
            return
        self._decrement_counter(self.read_ps_counts, old_value)
        self.read_ps_counts[new_value] += 1
        self.read_phase_set_by_index[index] = new_value
        self.read_arrays["ps"][index] = new_value

    def _register_reads_bulk(self, records: Iterable[dict[str, Any]], timing: dict[str, Any] | None = None) -> None:
        step_started = time.perf_counter()
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
                    self._replace_read_hp(existing, 0)
                    self.haplotype_conflict_count += 1
                old_phase_set = self.read_phase_set_by_index[existing]
                if old_phase_set != phase_set and phase_set != 0:
                    self._replace_read_ps(existing, 0)
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

        _timing_add(timing, "register_reads_dedup_conflict", time.perf_counter() - step_started)
        _timing_set_count(timing, "new_shared_reads", len(pending_order))
        if not pending_order:
            return

        old = len(self.read_lengths_by_index)
        sequences = [str(pending[key]["read_sequence"]).upper() for key in pending_order]
        step_started = time.perf_counter()
        packed, local_offsets, lengths, total_bases, invalid_counts = pack_4bit_acgtnp_chunk(
            sequences,
            READ_ALLOWED,
            old,
            "fail",
        )
        _timing_add(timing, "pack_read_sequences_4bit", time.perf_counter() - step_started)
        _timing_set_count(timing, "read_bases_packed", int(total_bases))
        if invalid_counts:
            raise ValueError(f"unexpected read invalid base counts: {invalid_counts}")
        if np.any(lengths <= 0):
            raise ValueError("empty shared read sequence")
        step_started = time.perf_counter()
        haplotypes = np.asarray([int(pending[key]["read_haplotype"]) for key in pending_order], dtype=np.uint8)
        phase_sets = np.asarray([int(pending[key].get("read_phase_set", 0)) for key in pending_order], dtype=np.int64)
        reference_starts = np.asarray([int(pending[key].get("reference_start", -1)) for key in pending_order], dtype=np.int64)
        reference_ends = np.asarray([int(pending[key].get("reference_end", -1)) for key in pending_order], dtype=np.int64)
        segment_types = np.asarray([int(pending[key].get("segment_type", 0)) for key in pending_order], dtype=np.int8)
        cigar_lists = [np.asarray(pending[key].get("cigar_ops_packed", []), dtype=np.uint32) for key in pending_order]
        cigar_lens = np.asarray([int(values.shape[0]) for values in cigar_lists], dtype=np.int32)
        cigar_offsets = self.cigar_op_count + np.concatenate(([0], np.cumsum(cigar_lens[:-1], dtype=np.int64))).astype(np.int64)
        if np.any(~np.isin(haplotypes, np.asarray([0, 1, 2], dtype=np.uint8))):
            raise ValueError("invalid HP tag value outside 0/1/2")
        if np.any(phase_sets < 0):
            raise ValueError("negative phase-set id")
        if np.any(reference_ends < reference_starts):
            raise ValueError("invalid reference interval")
        if np.any(segment_types < 0):
            raise ValueError("negative segment type")
        _timing_add(timing, "build_read_metadata_arrays", time.perf_counter() - step_started)

        step_started = time.perf_counter()
        old_base_offset = self.read_base_offset
        new_read_offsets = local_offsets[1:] + old_base_offset
        append_1d(self.read_arrays["bases_packed"], packed)
        append_1d(self.read_arrays["lengths"], lengths.astype(np.int32, copy=False))
        append_1d(self.read_arrays["offsets"], new_read_offsets)
        self.read_base_offset = old_base_offset + int(packed.shape[0])
        _timing_add_count(timing, "read_sequence_append_calls", 3)
        _timing_add(timing, "write_read_sequence_arrays", time.perf_counter() - step_started)

        step_started = time.perf_counter()
        total_new_cigar_ops = int(cigar_lens.sum())
        if total_new_cigar_ops:
            old_cigar = self.cigar_op_count
            new_cigar = old_cigar + total_new_cigar_ops
            self.cigar_ops.resize((new_cigar,))
            self.cigar_ops[old_cigar:new_cigar] = np.concatenate(cigar_lists).astype(np.uint32, copy=False)
            self.cigar_op_count = new_cigar
            _timing_add_count(timing, "cigar_array_writes", 1)
        _timing_add(timing, "write_cigar_arrays", time.perf_counter() - step_started)

        step_started = time.perf_counter()
        read_metadata_items = (
            ("hp", haplotypes),
            ("ps", phase_sets),
            ("ref_start", reference_starts),
            ("ref_end", reference_ends),
            ("segment_type", segment_types),
            ("cigar_offset", cigar_offsets),
            ("cigar_len", cigar_lens),
        )
        for name, values in read_metadata_items:
            append_1d(self.read_arrays[name], values)
        _timing_add_count(timing, "read_metadata_append_calls", len(read_metadata_items))
        _timing_add(timing, "write_read_metadata_arrays", time.perf_counter() - step_started)

        step_started = time.perf_counter()
        for offset, key in enumerate(pending_order, start=old):
            local_index = offset - old
            self.read_key_to_index[key] = offset
            self.read_haplotype_by_index.append(int(haplotypes[local_index]))
            self.read_phase_set_by_index.append(int(phase_sets[local_index]))
            self.read_lengths_by_index.append(int(lengths[local_index]))
            self.read_index.reference_ends.append(int(reference_ends[local_index]))
        self.read_hp_counts.update(int(value) for value in haplotypes)
        self.read_ps_counts.update(int(value) for value in phase_sets)
        batch_min = int(lengths.min())
        batch_max = int(lengths.max())
        self.read_length_min = batch_min if self.read_length_min is None else min(self.read_length_min, batch_min)
        self.read_length_max = batch_max if self.read_length_max is None else max(self.read_length_max, batch_max)
        self.read_length_sum += int(total_bases)
        _timing_add(timing, "update_read_indices", time.perf_counter() - step_started)

    def append_hap_chunk(self, label: str, artifacts: dict[str, Any], timing: dict[str, Any] | None = None) -> None:
        if label not in self.hap_arrays:
            raise ValueError(f"unknown hap label: {label}")
        windows = artifacts["windows"]
        records = artifacts["records"]
        if len(windows) != len(records):
            raise ValueError("windows/records length mismatch")
        if self.hap_window_cursor[label] + len(windows) > self.hap_window_counts[label]:
            raise ValueError(f"too many windows appended for {label}")

        step_started = time.perf_counter()
        input_segment_count = sum(len(window_records) for window_records in records)
        _timing_add(timing, "count_input_segments", time.perf_counter() - step_started)
        _timing_set_count(timing, "input_segments", input_segment_count)
        self._register_reads_bulk((record for window_records in records for record in window_records), timing=timing)
        step_started = time.perf_counter()
        self._write_hap_provenance_attrs(label, artifacts.get("attrs", {}))
        _timing_add(timing, "write_hap_provenance_attrs", time.perf_counter() - step_started)

        arrays = self.hap_arrays[label]
        wlo = self.hap_window_cursor[label]
        whi = wlo + len(windows)
        sequences = [str(w.sequence).upper() for w in windows]
        step_started = time.perf_counter()
        packed_seq, local_offsets, lengths, _total_bases, invalid_counts = pack_4bit_acgtnp_chunk(
            sequences,
            SEQUENCE_ALLOWED,
            wlo,
            "to_n",
        )
        _timing_add(timing, "pack_hap_sequences_4bit", time.perf_counter() - step_started)

        for window, sequence in zip(windows, sequences, strict=True):
            normalized = sequence if set(sequence) <= SEQUENCE_ALLOWED else "".join(
                base if base in SEQUENCE_ALLOWED else "N" for base in sequence
            )
            self.hap_window_content_hash[label].update(
                f"{int(window.raw_start)}:{int(window.actual_length)}:{len(normalized)}\\n".encode("ascii")
            )
            self.hap_window_content_hash[label].update(normalized.encode("ascii"))
        step_started = time.perf_counter()
        old_base_offset = self.hap_sequence_offsets[label][-1]
        new_hap_offsets = local_offsets[1:] + old_base_offset
        append_1d(arrays["bases_packed"], packed_seq)
        append_1d(arrays["lengths"], lengths.astype(np.int32, copy=False))
        append_1d(arrays["offsets"], new_hap_offsets)
        raw_starts = np.asarray([w.raw_start for w in windows], dtype=np.int64)
        actual_lengths = np.asarray([w.actual_length for w in windows], dtype=np.int64)
        arrays["raw_start"][wlo:whi] = raw_starts
        arrays["actual_length"][wlo:whi] = actual_lengths
        self.hap_sequence_offsets[label].extend(int(value) for value in new_hap_offsets)
        self.hap_raw_starts[label].extend(int(value) for value in raw_starts)
        self.hap_actual_lengths[label].extend(int(value) for value in actual_lengths)
        self.hap_total_bases[label] += int(lengths.sum())
        self.hap_invalid_counts[label].update(invalid_counts)
        _timing_add_count(timing, "hap_sequence_append_calls", 3)
        _timing_add_count(timing, "hap_window_metadata_writes", 2)
        _timing_add(timing, "write_hap_arrays", time.perf_counter() - step_started)

        step_started = time.perf_counter()
        depth_flat, depth_local_offsets, depth_lengths, depth_sum, depth_max = _depth_arrays_for_windows(
            windows,
            records,
            method=self.depth_method,
            segment_types=self.depth_segment_types,
        )
        if depth_lengths.shape[0] != lengths.shape[0] or np.any(depth_lengths != lengths.astype(np.int32, copy=False)):
            raise ValueError(f"{label} depth/window length mismatch")
        if not hasattr(self, "window_quality"):
            self.window_quality = []
        for i, window in enumerate(windows):
            begin = int(depth_local_offsets[i])
            actual = int(window.actual_length)
            self.window_quality.append({"contig": self.contig, "haplotype": label,
                "window_index": wlo+i, "raw_start": int(window.raw_start), "actual_length": actual,
                **depth_metrics(depth_flat[begin:begin+actual], actual)})
        _timing_add(timing, "build_depth_arrays", time.perf_counter() - step_started)

        step_started = time.perf_counter()
        depth_arrays = self.depth_arrays[label]
        old_depth_offset = self.depth_offsets[label][-1]
        new_depth_offsets = depth_local_offsets[1:] + old_depth_offset
        append_1d(depth_arrays["depth"], depth_flat)
        append_1d(depth_arrays["lengths"], depth_lengths)
        append_1d(depth_arrays["offsets"], new_depth_offsets)
        self.depth_offsets[label].extend(int(value) for value in new_depth_offsets)
        self.depth_total_bases[label] += int(depth_lengths.sum())
        self.depth_sum[label] += int(depth_sum)
        self.depth_max[label] = max(self.depth_max[label], int(depth_max))
        _timing_add_count(timing, "depth_append_calls", 3)
        _timing_add(timing, "write_depth_arrays", time.perf_counter() - step_started)
        if self.has_supplementary_depth:
            step_started = time.perf_counter()
            supplementary_flat, supplementary_local_offsets, supplementary_lengths, supplementary_sum, supplementary_max = _depth_arrays_for_windows(
                windows, records, method=self.depth_method, segment_types=self.supplementary_depth_segment_types,
            )
            if supplementary_lengths.shape[0] != lengths.shape[0] or np.any(supplementary_lengths != lengths.astype(np.int32, copy=False)):
                raise ValueError(f"{label} supplementary depth/window length mismatch")
            supplementary_arrays = self.supplementary_depth_arrays[label]
            old_supplementary_offset = self.supplementary_depth_offsets[label][-1]
            new_supplementary_offsets = supplementary_local_offsets[1:] + old_supplementary_offset
            append_1d(supplementary_arrays["depth"], supplementary_flat)
            append_1d(supplementary_arrays["lengths"], supplementary_lengths)
            append_1d(supplementary_arrays["offsets"], new_supplementary_offsets)
            self.supplementary_depth_offsets[label].extend(int(value) for value in new_supplementary_offsets)
            self.supplementary_depth_total_bases[label] += int(supplementary_lengths.sum())
            self.supplementary_depth_sum[label] += int(supplementary_sum)
            self.supplementary_depth_max[label] = max(self.supplementary_depth_max[label], int(supplementary_max))
            _timing_add_count(timing, "supplementary_depth_append_calls", 3)
            _timing_add(timing, "build_write_supplementary_depth_arrays", time.perf_counter() - step_started)


        step_started = time.perf_counter()
        packed_segments: list[tuple[int, int, int, int, int, int]] = []
        local_segment_offsets = [0]
        for window, window_records in zip(windows, records):
            row_values: list[tuple[int, int, int, int, int, int]] = []
            for record in window_records:
                key = tuple(record.get("shared_read_key", record["read_key"]))
                read_index = self.read_key_to_index[key]
                query_start = int(record["query_start"])
                query_end = int(record["query_end"])
                hap_start = int(record["read_start"])
                hap_end = int(record["read_end"])
                read_strand = int(record["read_strand"])
                if query_start < 0 or query_end <= query_start:
                    raise ValueError("invalid query interval")
                if read_index < 0 or read_index >= len(self.read_lengths_by_index):
                    raise ValueError("shared read index out of range")
                read_length = self.read_lengths_by_index[read_index]
                if query_end > read_length:
                    raise ValueError(f"query interval outside shared read: {query_start}-{query_end} length={read_length}")
                if hap_start < 0 or hap_end <= hap_start or hap_end > int(window.actual_length):
                    raise ValueError(f"invalid hap read interval: {hap_start}-{hap_end}")
                if read_strand not in (-1, 1):
                    raise ValueError("invalid read strand")
                row_values.append((read_index, query_start, query_end, hap_start, hap_end, read_strand))
            row_values.sort(key=lambda x: (x[0], x[1], x[2], x[3], x[4]))
            packed_segments.extend(row_values)
            local_segment_offsets.append(len(packed_segments))

        _timing_add(timing, "build_segment_arrays", time.perf_counter() - step_started)
        _timing_set_count(timing, "segments_written", len(packed_segments))
        step_started = time.perf_counter()
        old_segment = self.hap_segment_cursor[label]
        new_segment = old_segment + len(packed_segments)
        seg_arrays = self.segment_arrays[label]
        if packed_segments:
            values = np.asarray(packed_segments, dtype=np.int64)
            append_1d(seg_arrays["read_id"], values[:, 0].astype(np.int32, copy=False))
            append_1d(seg_arrays["query_start"], values[:, 1].astype(np.int32, copy=False))
            append_1d(seg_arrays["query_end"], values[:, 2].astype(np.int32, copy=False))
            append_1d(seg_arrays["hap_start"], values[:, 3].astype(np.int64, copy=False))
            append_1d(seg_arrays["hap_end"], values[:, 4].astype(np.int64, copy=False))
            append_1d(seg_arrays["strand"], values[:, 5].astype(np.int8, copy=False))
            _timing_add_count(timing, "segment_append_calls", 6)
        new_window_offsets = np.asarray(local_segment_offsets[1:], dtype=np.int64) + old_segment
        seg_arrays["window_read_offsets"][wlo + 1:whi + 1] = new_window_offsets
        _timing_add_count(timing, "window_read_offset_writes", 1 if len(windows) else 0)
        self.hap_segment_offsets[label].extend(int(value) for value in new_window_offsets)

        _timing_add(timing, "write_segment_arrays", time.perf_counter() - step_started)

        self.hap_window_cursor[label] = whi
        self.hap_segment_cursor[label] = new_segment
        self.hap_counts[label].update(artifacts.get("counts", {}))
        self.hap_coverage[label].extend(float(x) for x in artifacts.get("attrs", {}).get("selected_coordinate_coverage", []))

    def _write_hap_provenance_attrs(self, label: str, attrs: dict[str, Any]) -> None:
        if label in self.hap_provenance_written:
            return
        group = self.hap_groups[label]
        for key in HAP_PROVENANCE_ATTR_KEYS:
            if key in attrs and key not in group.attrs:
                group.attrs[key] = attrs[key]
        self.hap_provenance_written.add(label)

    def finalize(self, *, validation_level: str, timing: DirectTiming | None = None) -> dict[str, Any]:
        for label, total in self.hap_window_counts.items():
            if self.hap_window_cursor[label] != total:
                raise ValueError(f"window count mismatch for {label}: {self.hap_window_cursor[label]} != {total}")
        validation_summary: dict[str, Any] = {
            "schema": PERSONAL_FULL_VALIDATION_SCHEMA,
            "level": validation_level,
            "status": "pass",
            "checks": {"structural": True},
        }
        step_started = time.perf_counter()
        self._validate_shapes(validation_level)
        if validation_level == "full":
            validation_summary["checks"]["semantic"] = self._validate_full_semantics()
        if timing is not None:
            timing.add("finalize_validate_shapes", time.perf_counter() - step_started)
        if any(os.environ.get('MFM_'+key+'_BACKEND','cpu')=='cupy' for key in ['COVERAGE','QC','ZARR']):
            from .gpu_arrays import execution_metadata
            execution=execution_metadata()
            for group in [self.hap_seq_root,self.depth_root,self.readset_root]:group.attrs['compute_execution']=execution
        step_started = time.perf_counter()
        hap_summaries = self._hap_summaries()
        read_summary = self._read_summary()
        depth_summary = self._depth_summary()
        if timing is not None:
            timing.add("finalize_build_summaries", time.perf_counter() - step_started)
        step_started = time.perf_counter()
        windows = []
        for label, summary in hap_summaries.items():
            windows.extend(make_window_rows(self.contig, label, summary))
        write_windows_tsv(self.target_root / "windows.tsv", windows)
        write_quality(self.target_root / "window_quality.tsv", getattr(self, "window_quality", []))
        if timing is not None:
            timing.add("finalize_write_windows_tsv", time.perf_counter() - step_started)
        payload = {
            "schema": FULL_DIRECT_CONTIG_SCHEMA,
            "status": "direct_compact",
            "sample_id": self.sample,
            "source_sample_id": self.sample_metadata.get("source_sample_id", self.sample),
            "dataset_role": self.sample_metadata.get("dataset_role", ""),
            "read_type": self.sample_metadata.get("read_type", ""),
            "project_id": self.sample_metadata.get("project_id", ""),
            "layout_name": self.sample_metadata.get("layout_name", ""),
            "contig": self.contig,
            "source_alignment": str(self.source_alignment),
            "reference_fasta": str(self.reference_fasta) if self.reference_fasta else "",
            "target": str(self.final_target),
            "window_rows": len(windows),
            "hap_summaries": hap_summaries,
            "read_summary": read_summary,
            "depth_summary": depth_summary,
            "filter_policy": {
                "preset": self.filter_preset,
                "depth_segment_types": list(self.depth_segment_types),
                "supplementary_depth_segment_types": list(self.supplementary_depth_segment_types),
                "model_segment_types": list(self.model_segment_types),
            },
            "filter_qc": self.filter_qc,
            "encoding": "4bit_acgtnp",
            "hap_invalid_policy": "to_n",
            "validation": validation_summary,
            "hap_invalid_replacement": "N",
            "validation_level": validation_level,
            "completed_at": utc_now(),
        }
        if timing is not None:
            payload["direct_timing"] = timing.to_dict()
        step_started = time.perf_counter()
        write_json(self.target_root / "contig_metadata.json", payload)
        write_json(self.target_root / "filter_qc.json", self.filter_qc)
        if timing is not None:
            timing.add("finalize_write_contig_metadata", time.perf_counter() - step_started)
        return payload

    def _validate_shapes(self, validation_level: str) -> None:
        if validation_level not in {"full", "standard", "fast"}:
            raise ValueError(f"unknown validation level: {validation_level}")
        n_reads = int(self.read_arrays["lengths"].shape[0])
        if len(self.read_key_to_index) > n_reads or n_reads != len(self.read_lengths_by_index):
            raise ValueError("invalid in-memory read index length")
        if int(self.read_arrays["offsets"].shape[0]) != n_reads + 1:
            raise ValueError("invalid read offsets length")
        if int(self.read_arrays["hp"].shape[0]) != n_reads or int(self.read_arrays["ps"].shape[0]) != n_reads:
            raise ValueError("invalid read HP/PS lengths")
        for name in ("ref_start", "ref_end", "segment_type", "cigar_offset", "cigar_len"):
            if int(self.read_arrays[name].shape[0]) != n_reads:
                raise ValueError(f"invalid read {name} length")
        if int(self.read_arrays["bases_packed"].shape[0]) != int(self.read_base_offset):
            raise ValueError("invalid read packed base cursor")
        if int(self.cigar_ops.shape[0]) != int(self.cigar_op_count):
            raise ValueError("invalid CIGAR op cursor")

        for label in self.hap_window_counts:
            hap_arrays = self.hap_arrays[label]
            seg_arrays = self.segment_arrays[label]
            depth_arrays = self.depth_arrays[label]
            n_windows = int(self.hap_window_counts[label])
            n_segments = int(seg_arrays["read_id"].shape[0])
            if int(hap_arrays["offsets"].shape[0]) != n_windows + 1:
                raise ValueError(f"invalid {label} hap offsets length")
            if int(hap_arrays["lengths"].shape[0]) != n_windows:
                raise ValueError(f"invalid {label} hap lengths")
            if int(hap_arrays["raw_start"].shape[0]) != n_windows or int(hap_arrays["actual_length"].shape[0]) != n_windows:
                raise ValueError(f"invalid {label} hap window metadata length")
            if len(self.hap_sequence_offsets[label]) != n_windows + 1:
                raise ValueError(f"invalid {label} in-memory hap offsets length")
            if len(self.hap_raw_starts[label]) != n_windows or len(self.hap_actual_lengths[label]) != n_windows:
                raise ValueError(f"invalid {label} in-memory hap window metadata length")
            if int(hap_arrays["bases_packed"].shape[0]) != int(self.hap_sequence_offsets[label][-1]):
                raise ValueError(f"invalid {label} hap packed base cursor")
            if int(depth_arrays["offsets"].shape[0]) != n_windows + 1:
                raise ValueError(f"invalid {label} depth offsets length")
            if int(depth_arrays["lengths"].shape[0]) != n_windows:
                raise ValueError(f"invalid {label} depth lengths")
            if len(self.depth_offsets[label]) != n_windows + 1:
                raise ValueError(f"invalid {label} in-memory depth offsets length")
            if int(depth_arrays["depth"].shape[0]) != int(self.depth_offsets[label][-1]):
                raise ValueError(f"invalid {label} depth cursor")
            if self.has_supplementary_depth:
                supplementary_arrays = self.supplementary_depth_arrays[label]
                if int(supplementary_arrays["offsets"].shape[0]) != n_windows + 1:
                    raise ValueError(f"invalid {label} supplementary depth offsets length")
                if int(supplementary_arrays["lengths"].shape[0]) != n_windows:
                    raise ValueError(f"invalid {label} supplementary depth lengths")
                if len(self.supplementary_depth_offsets[label]) != n_windows + 1:
                    raise ValueError(f"invalid {label} in-memory supplementary depth offsets length")
                if int(supplementary_arrays["depth"].shape[0]) != int(self.supplementary_depth_offsets[label][-1]):
                    raise ValueError(f"invalid {label} supplementary depth cursor")
            if int(seg_arrays["window_read_offsets"].shape[0]) != n_windows + 1:
                raise ValueError(f"invalid {label} window segment offsets length")
            if len(self.hap_segment_offsets[label]) != n_windows + 1:
                raise ValueError(f"invalid {label} in-memory window segment offsets length")
            if int(self.hap_segment_offsets[label][-1]) != n_segments:
                raise ValueError(f"invalid {label} segment cursor")
            for name in ("hap_start", "hap_end", "query_start", "query_end", "strand"):
                if int(seg_arrays[name].shape[0]) != n_segments:
                    raise ValueError(f"invalid {label} segment {name} length")

        if validation_level == "fast":
            return

        offsets = self.read_arrays["offsets"][:]
        if offsets.size and (int(offsets[0]) != 0 or int(offsets[-1]) != int(self.read_base_offset) or np.any(np.diff(offsets) < 0)):
            raise ValueError("invalid read offsets")
        lengths = self.read_arrays["lengths"][:]
        if n_reads and np.any(lengths <= 0):
            raise ValueError("invalid read lengths")
        if np.any(~np.isin(self.read_arrays["hp"][:], np.asarray([0, 1, 2], dtype=np.uint8))):
            raise ValueError("invalid HP values")
        if np.any(self.read_arrays["ps"][:] < 0):
            raise ValueError("invalid PS values")
        if np.any(self.read_arrays["ref_end"][:] < self.read_arrays["ref_start"][:]):
            raise ValueError("invalid reference intervals")
        cigar_offsets = self.read_arrays["cigar_offset"][:]
        cigar_lens = self.read_arrays["cigar_len"][:]
        if np.any(cigar_offsets < 0) or np.any(cigar_lens < 0) or np.any(cigar_offsets + cigar_lens > int(self.cigar_ops.shape[0])):
            raise ValueError("invalid CIGAR offsets")
        for label in self.hap_window_counts:
            hap_arrays = self.hap_arrays[label]
            seg_arrays = self.segment_arrays[label]
            depth_arrays = self.depth_arrays[label]
            n_segments = int(seg_arrays["read_id"].shape[0])
            depth_offsets = depth_arrays["offsets"][:]
            depth_lengths = depth_arrays["lengths"][:]
            if int(depth_offsets[0]) != 0 or int(depth_offsets[-1]) != int(self.depth_offsets[label][-1]) or np.any(np.diff(depth_offsets) < 0):
                raise ValueError(f"invalid {label} depth offsets")
            if np.any(depth_lengths < 0):
                raise ValueError(f"invalid {label} depth lengths")
            if depth_lengths.shape[0] and np.any(depth_lengths != hap_arrays["lengths"][:]):
                raise ValueError(f"invalid {label} depth lengths do not match hap lengths")
            if np.any(depth_offsets[1:] - depth_offsets[:-1] != depth_lengths):
                raise ValueError(f"invalid {label} depth offset deltas")
            if self.has_supplementary_depth:
                supplementary_arrays = self.supplementary_depth_arrays[label]
                supplementary_offsets = supplementary_arrays["offsets"][:]
                supplementary_lengths = supplementary_arrays["lengths"][:]
                if int(supplementary_offsets[0]) != 0 or int(supplementary_offsets[-1]) != int(self.supplementary_depth_offsets[label][-1]) or np.any(np.diff(supplementary_offsets) < 0):
                    raise ValueError(f"invalid {label} supplementary depth offsets")
                if np.any(supplementary_lengths < 0):
                    raise ValueError(f"invalid {label} supplementary depth lengths")
                if supplementary_lengths.shape[0] and np.any(supplementary_lengths != hap_arrays["lengths"][:]):
                    raise ValueError(f"invalid {label} supplementary depth lengths do not match hap lengths")
                if np.any(supplementary_offsets[1:] - supplementary_offsets[:-1] != supplementary_lengths):
                    raise ValueError(f"invalid {label} supplementary depth offset deltas")
            window_offsets = seg_arrays["window_read_offsets"][:]
            if int(window_offsets[0]) != 0 or int(window_offsets[-1]) != n_segments or np.any(np.diff(window_offsets) < 0):
                raise ValueError(f"invalid {label} window segment offsets")
            hap_offsets = hap_arrays["offsets"][:]
            if int(hap_offsets[0]) != 0 or int(hap_offsets[-1]) != int(self.hap_sequence_offsets[label][-1]) or np.any(np.diff(hap_offsets) < 0):
                raise ValueError(f"invalid {label} hap offsets")
            if n_segments:
                if np.any(seg_arrays["read_id"][:] < 0) or np.any(seg_arrays["read_id"][:] >= n_reads):
                    raise ValueError(f"invalid {label} read ids")
                if np.any(seg_arrays["query_end"][:] <= seg_arrays["query_start"][:]):
                    raise ValueError(f"invalid {label} query intervals")
                if np.any(seg_arrays["hap_end"][:] <= seg_arrays["hap_start"][:]):
                    raise ValueError(f"invalid {label} hap intervals")
                if np.any(~np.isin(seg_arrays["strand"][:], np.asarray([-1, 1], dtype=np.int8))):
                    raise ValueError(f"invalid {label} strand values")


    def _validate_packed_4bit(self, *, name: str, bases, offsets, lengths, max_code: int) -> dict[str, int]:
        record_count = int(lengths.shape[0])
        byte_count = int(bases.shape[0])
        for lo in range(0, record_count, VALIDATION_RECORD_CHUNK):
            hi = min(record_count, lo + VALIDATION_RECORD_CHUNK)
            block_offsets = np.asarray(offsets[lo:hi + 1], dtype=np.int64)
            block_lengths = np.asarray(lengths[lo:hi], dtype=np.int64)
            expected = (block_lengths + 1) // 2
            if np.any(block_lengths <= 0) or np.any(np.diff(block_offsets) != expected):
                raise ValueError(f"full validation invalid packed {name} offsets/lengths at records {lo}:{hi}")
        for lo in range(0, byte_count, VALIDATION_PACKED_BYTE_CHUNK):
            hi = min(byte_count, lo + VALIDATION_PACKED_BYTE_CHUNK)
            values = np.asarray(bases[lo:hi], dtype=np.uint8)
            if np.any((values >> 4) > max_code) or np.any((values & 0x0F) > max_code):
                raise ValueError(f"full validation invalid 4-bit {name} code at bytes {lo}:{hi}")
        return {"records": record_count, "packed_bytes": byte_count}

    def _validate_full_semantics(self) -> dict[str, Any]:
        """Recalculate all encoded CIGAR/depth relationships without source BAM reads."""
        n_reads = int(self.read_arrays["lengths"].shape[0])
        read_cache = _ReadValidationChunkCache(
            self.read_arrays, self.cigar_ops, n_reads, chunk_size=VALIDATION_RECORD_CHUNK
        )
        packed_read = self._validate_packed_4bit(
            name="read",
            bases=self.read_arrays["bases_packed"],
            offsets=self.read_arrays["offsets"],
            lengths=self.read_arrays["lengths"],
            max_code=4,
        )
        cigar_backend = os.environ.get('MFM_QC_CIGAR_BACKEND','numpy')
        if cigar_backend not in {'numpy','legacy'}:
            raise ValueError('MFM_QC_CIGAR_BACKEND는 numpy|legacy여야 함')
        if cigar_backend == 'legacy':
            for read_id in range(n_reads):read_cache.decode(read_id)
            cigar_execution={'mode':'legacy','records':n_reads,'threads':1}
        else:
            from .compact.chunk_validation import validate_cigars
            cigar_execution=validate_cigars(self.read_arrays,self.cigar_ops,n_reads)
        print('[QC-CIGAR] '+json.dumps(cigar_execution),flush=True)

        depth_windows = 0
        depth_bases = 0
        supplementary_depth_windows = 0
        packed_hap: dict[str, dict[str, Any]] = {}
        for label in self.hap_window_counts:
            hap_arrays = self.hap_arrays[label]
            seg_arrays = self.segment_arrays[label]
            depth_arrays = self.depth_arrays[label]
            supplementary_arrays = self.supplementary_depth_arrays.get(label)
            supplementary_offsets = np.asarray(supplementary_arrays["offsets"][:], dtype=np.int64) if supplementary_arrays is not None else None
            packed_hap[label] = self._validate_packed_4bit(
                name=f"{label} hap",
                bases=hap_arrays["bases_packed"],
                offsets=hap_arrays["offsets"],
                lengths=hap_arrays["lengths"],
                max_code=5,
            )
            raw_starts = np.asarray(hap_arrays["raw_start"][:], dtype=np.int64)
            actual_lengths = np.asarray(hap_arrays["actual_length"][:], dtype=np.int64)
            storage_lengths = np.asarray(hap_arrays["lengths"][:], dtype=np.int64)
            window_offsets = np.asarray(seg_arrays["window_read_offsets"][:], dtype=np.int64)
            depth_offsets = np.asarray(depth_arrays["offsets"][:], dtype=np.int64)
            hap_offsets = np.asarray(hap_arrays["offsets"][:], dtype=np.int64)
            observed_window_hash = hashlib.sha256()
            for window_index, (raw_start, actual_length, storage_length) in enumerate(
                zip(raw_starts, actual_lengths, storage_lengths, strict=True)
            ):
                if actual_length <= 0 or actual_length > storage_length:
                    raise ValueError(f"full validation invalid {label} window length at index={window_index}")
                byte_lo = int(hap_offsets[window_index])
                byte_hi = int(hap_offsets[window_index + 1])
                packed_window = np.asarray(hap_arrays["bases_packed"][byte_lo:byte_hi], dtype=np.uint8)
                codes = np.empty(packed_window.size * 2, dtype=np.uint8)
                codes[0::2] = packed_window >> 4
                codes[1::2] = packed_window & 0x0F
                sequence = FOURBIT_CODE_TO_ASCII[codes[:int(storage_length)]].tobytes()
                actual_sequence = sequence[:int(actual_length)]
                if set(actual_sequence) - set(b"ACGTN") or sequence[int(actual_length):] != b"P" * (int(storage_length) - int(actual_length)):
                    raise ValueError(f"invalid {label} HAP N/padding policy at window={window_index}")
                if not any(base in b"ACGT" for base in actual_sequence):
                    raise ValueError(f"uninformative all-N {label} window={window_index}")
                observed_window_hash.update(
                    f"{int(raw_start)}:{int(actual_length)}:{int(storage_length)}\\n".encode("ascii")
                )
                observed_window_hash.update(sequence)
                segment_lo = int(window_offsets[window_index])
                segment_hi = int(window_offsets[window_index + 1])
                read_ids = np.asarray(seg_arrays["read_id"][segment_lo:segment_hi], dtype=np.int64)
                query_starts = np.asarray(seg_arrays["query_start"][segment_lo:segment_hi], dtype=np.int64)
                query_ends = np.asarray(seg_arrays["query_end"][segment_lo:segment_hi], dtype=np.int64)
                hap_starts = np.asarray(seg_arrays["hap_start"][segment_lo:segment_hi], dtype=np.int64)
                hap_ends = np.asarray(seg_arrays["hap_end"][segment_lo:segment_hi], dtype=np.int64)
                if read_ids.size and np.any(np.diff(read_ids) <= 0):
                    raise ValueError(f"full validation duplicate/unsorted {label} read ids at window={window_index}")
                diff = np.zeros(int(storage_length) + 1, dtype=np.int64)
                supplementary_diff = np.zeros(int(storage_length) + 1, dtype=np.int64) if supplementary_arrays is not None else None
                window_end = int(raw_start + actual_length)
                for read_id, query_start, query_end, hap_start, hap_end in zip(
                    read_ids, query_starts, query_ends, hap_starts, hap_ends, strict=True
                ):
                    read_index = int(read_id)
                    intervals, _query, _reference, read_length, segment_type = read_cache.decode(read_index)
                    if query_start < 0 or query_end <= query_start or query_end > read_length:
                        raise ValueError(f"full validation invalid {label} query segment at window={window_index}")
                    if hap_start < 0 or hap_end <= hap_start or hap_end > actual_length:
                        raise ValueError(f"full validation invalid {label} HAP segment at window={window_index}")
                    target_diffs = []
                    if segment_type in self.depth_segment_types:
                        target_diffs.append(diff)
                    if supplementary_diff is not None and segment_type in self.supplementary_depth_segment_types:
                        target_diffs.append(supplementary_diff)
                    for ref_start, ref_end in intervals:
                        start = max(int(ref_start), int(raw_start))
                        end = min(int(ref_end), window_end)
                        if end > start:
                            for target_diff in target_diffs:
                                target_diff[start - int(raw_start)] += 1
                                target_diff[end - int(raw_start)] -= 1
                expected_depth = np.cumsum(diff[:-1], dtype=np.int64).astype(np.uint32, copy=False)
                depth_lo = int(depth_offsets[window_index])
                depth_hi = int(depth_offsets[window_index + 1])
                observed_depth = np.asarray(depth_arrays["depth"][depth_lo:depth_hi], dtype=np.uint32)
                if observed_depth.shape[0] != int(storage_length) or not np.array_equal(observed_depth, expected_depth):
                    raise ValueError(f"full validation CIGAR-derived depth mismatch: {label} window={window_index}")
                if supplementary_arrays is not None and supplementary_offsets is not None and supplementary_diff is not None:
                    expected_supplementary = np.cumsum(supplementary_diff[:-1], dtype=np.int64).astype(np.uint32, copy=False)
                    supplementary_lo = int(supplementary_offsets[window_index])
                    supplementary_hi = int(supplementary_offsets[window_index + 1])
                    observed_supplementary = np.asarray(supplementary_arrays["depth"][supplementary_lo:supplementary_hi], dtype=np.uint32)
                    if observed_supplementary.shape[0] != int(storage_length) or not np.array_equal(observed_supplementary, expected_supplementary):
                        raise ValueError(f"full validation CIGAR-derived supplementary depth mismatch: {label} window={window_index}")
                    supplementary_depth_windows += 1

                depth_windows += 1
                depth_bases += int(storage_length)

            observed_digest = observed_window_hash.hexdigest()
            if observed_digest != self.hap_window_content_hash[label].hexdigest():
                raise ValueError(f"full validation HAP window sequence mismatch: {label}")
            packed_hap[label]["window_content_sha256"] = observed_digest

        return {
            "read_count": n_reads,
            "cigar_op_count": int(self.cigar_ops.shape[0]),
            "packed_read": packed_read,
            "packed_haplotype": packed_hap,
            "depth_windows": depth_windows,
            "depth_bases": depth_bases,
            "supplementary_depth_windows": supplementary_depth_windows,
        }


    def _hap_summaries(self) -> dict[str, Any]:
        summaries = {}
        for label in self.hap_window_counts:
            offsets = self.hap_sequence_offsets[label]
            raw_starts = self.hap_raw_starts[label]
            actual_lengths = self.hap_actual_lengths[label]
            if len(offsets) != int(self.hap_window_counts[label]) + 1:
                raise ValueError(f"invalid {label} summary offsets length")
            summary = {
                "hap": label,
                "window_count": int(self.hap_window_counts[label]),
                "segment_count": int(self.hap_segment_cursor[label]),
                "total_bases": int(self.hap_total_bases[label]),
                "packed_bytes": int(offsets[-1]),
                "invalid_base_policy": "to_n",
                "invalid_base_count": int(sum(self.hap_invalid_counts[label].values())),
                "invalid_base_symbols": dict(sorted(self.hap_invalid_counts[label].items())),
                "raw_start": list(raw_starts),
                "actual_length": list(actual_lengths),
                "sequence_offsets": list(offsets),
                "window_content_sha256": self.hap_window_content_hash[label].hexdigest(),
                "filter_counts": dict(sorted(self.hap_counts[label].items())),
                "coordinate_coverage_min": min(self.hap_coverage[label]) if self.hap_coverage[label] else None,
                "coordinate_coverage_max": max(self.hap_coverage[label]) if self.hap_coverage[label] else None,
            }
            for key in HAP_PROVENANCE_ATTR_KEYS:
                if key in self.hap_groups[label].attrs:
                    summary[key] = self.hap_groups[label].attrs[key]
            summaries[label] = summary
        return summaries

    def _depth_summary(self) -> dict[str, Any]:
        by_haplotype: dict[str, Any] = {}
        total_bases = 0
        total_depth_sum = 0
        for label in self.hap_window_counts:
            bases = int(self.depth_total_bases[label])
            depth_sum = int(self.depth_sum[label])
            total_bases += bases
            total_depth_sum += depth_sum
            by_haplotype[label] = {
                "window_count": int(self.hap_window_counts[label]),
                "total_bases": bases,
                "depth_sum": depth_sum,
                "mean_depth": (depth_sum / bases) if bases else None,
                "max_depth": int(self.depth_max[label]),
                "depth_offsets": list(self.depth_offsets[label]),
                "depth_method": self.depth_method,
            }
        supplementary_evidence: dict[str, Any] | None = None
        if self.has_supplementary_depth:
            evidence_by_haplotype: dict[str, Any] = {}
            evidence_total_bases = 0
            evidence_depth_sum = 0
            for label in self.hap_window_counts:
                bases = int(self.supplementary_depth_total_bases[label])
                depth_sum = int(self.supplementary_depth_sum[label])
                evidence_total_bases += bases
                evidence_depth_sum += depth_sum
                evidence_by_haplotype[label] = {
                    "window_count": int(self.hap_window_counts[label]),
                    "total_bases": bases,
                    "depth_sum": depth_sum,
                    "mean_depth": (depth_sum / bases) if bases else None,
                    "max_depth": int(self.supplementary_depth_max[label]),
                    "depth_offsets": list(self.supplementary_depth_offsets[label]),
                    "depth_method": self.depth_method,
                }
            supplementary_evidence = {
                "layout": "per_haplotype_flat_uint32",
                "input_segment_types": list(self.supplementary_depth_segment_types),
                "total_bases": evidence_total_bases,
                "depth_sum": evidence_depth_sum,
                "mean_depth": (evidence_depth_sum / evidence_total_bases) if evidence_total_bases else None,
                "by_haplotype": evidence_by_haplotype,
            }

        return {
            "layout": "per_haplotype_flat_uint32",
            "dtype": "uint32",
            "method": self.depth_method,
            "read_coordinate_mode": self.read_coordinate_mode,
            "input_segment_types": list(self.depth_segment_types),
            "total_bases": total_bases,
            "depth_sum": total_depth_sum,
            "mean_depth": (total_depth_sum / total_bases) if total_bases else None,
            "by_haplotype": by_haplotype,
            "supplementary_evidence": supplementary_evidence,
        }

    def _read_summary(self) -> dict[str, Any]:
        n_reads = len(self.read_lengths_by_index)
        hp_counts = {str(int(value)): int(count) for value, count in sorted(self.read_hp_counts.items())}
        nonzero_ps_count = sum(int(count) for value, count in self.read_ps_counts.items() if int(value) > 0)
        zero_ps_count = int(self.read_ps_counts.get(0, 0))
        unique_nonzero_ps_count = sum(1 for value, count in self.read_ps_counts.items() if int(value) > 0 and int(count) > 0)
        segment_type_counts = Counter(int(value) for value in self.read_arrays["segment_type"][:])
        return {
            "read_count": n_reads,
            "total_bases": int(self.read_length_sum),
            "packed_bytes": int(self.read_base_offset),
            "invalid_base_count": 0,
            "invalid_base_symbols": {},
            "has_cigar": True,
            "cigar_op_count": int(self.cigar_op_count),
            "read_length_min": self.read_length_min,
            "read_length_max": self.read_length_max,
            "read_length_mean": (self.read_length_sum / n_reads) if n_reads else None,
            "hp_counts": hp_counts,
            "ps_summary": {
                "nonzero_count": int(nonzero_ps_count),
                "zero_count": int(zero_ps_count),
                "unique_nonzero_count": int(unique_nonzero_ps_count),
            },
            "segment_type_counts": {
                "primary": int(segment_type_counts.get(0, 0)),
                "supplementary": int(segment_type_counts.get(1, 0)),
            },
            "haplotype_conflict_count": int(self.haplotype_conflict_count),
            "phase_set_conflict_count": int(self.phase_set_conflict_count),
        }


def make_window_rows(contig: str, hap: str, summary: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for i in range(summary["window_count"]):
        start = int(summary["sequence_offsets"][i])
        end = int(summary["sequence_offsets"][i + 1])
        rows.append(
            {
                "contig": contig,
                "haplotype": hap,
                "window_index": i,
                "raw_start": int(summary["raw_start"][i]),
                "actual_length": int(summary["actual_length"][i]),
                "sequence_offset": start,
                "sequence_packed_bytes": end - start,
            }
        )
    return rows


def merge_count_dicts(items: list[dict[str, int]]) -> dict[str, int]:
    merged: dict[str, int] = {}
    for item in items:
        for key, value in item.items():
            merged[key] = merged.get(key, 0) + int(value)
    return dict(sorted(merged.items()))
