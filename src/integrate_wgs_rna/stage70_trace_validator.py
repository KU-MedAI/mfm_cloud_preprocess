"""Stage 70 personal HAP CRAM→Zarr 독립 추적 검증기.

writer/artifact generator를 호출하지 않고 source CRAM과 Zarr을 직접 decode한다.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import sqlite3
import tempfile
import time
from collections import Counter
from .cpu_executor import cpu_results, elastic_cpu_enabled
from multiprocessing import get_context
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pysam
import zarr
from .compact.runtime import positive_env_int, nonnegative_env_int
from .trace_acceleration import DepthAccumulator, WindowLookup
from .compact.trace_store import ReplayRecords
from .compact.trace_hash import digests,source_digests,hash_backend,segment_digests,SegmentBatch,add_hashes
from .compact.trace_decode import decode_reads,read_digest,segment_digest
from .compact.trace_multiset import PartitionedBag
from .compact.trace_resources import scratch_estimate,select_backend
from .window_chunk_cache import DecodedChunkCache
from .memory_policy import GIB, memory_limit

from numcodecs import blosc

HAPS = ("HAP1", "HAP2")
TRACE_SCHEMA = "stage70_source_to_zarr_trace_v1"
READ_BASES = np.frombuffer(b"ACGTN", dtype=np.uint8)
HAP_BASES = np.frombuffer(b"ACGTNP", dtype=np.uint8)
SQLITE_MULTISET_BATCH = positive_env_int("MFM_TRACE_SQLITE_BATCH", 8192)
ZARR_PREFETCH_RECORDS = positive_env_int("MFM_TRACE_PREFETCH_RECORDS", 65536)
SQLITE_CACHE_MIB = positive_env_int("MFM_TRACE_SQLITE_CACHE_MIB", 64)
# trace는 compact writer와 별도로 평면 contig×HAP pool을 쓴다. worker 수와
# Blosc thread 곱은 config_v2.sh의 sample CPU budget 아래로 제한된다.
TRACE_BLOSC_THREADS = positive_env_int(
    "MFM_TRACE_BLOSC_THREADS", positive_env_int("MFM_BLOSC_THREADS", 8)
)
blosc.set_nthreads(TRACE_BLOSC_THREADS)


def _json_hash(value: Any) -> bytes:
    return hashlib.sha256(
        json.dumps(value, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).digest()


def _unpack(values: Any, length: int, alphabet: np.ndarray, name: str) -> str:
    packed = np.asarray(values, dtype=np.uint8)
    codes = np.empty(packed.size * 2, dtype=np.uint8)
    codes[0::2], codes[1::2] = packed >> 4, packed & 15
    if length < 0 or length > codes.size or (codes[:length].size and codes[:length].max() >= alphabet.size):
        raise ValueError(f"invalid packed {name}")
    return alphabet[codes[:length]].tobytes().decode("ascii")


class _Bag:
    """Disk-backed multiset of per-record SHA-256 identities, never aggregate XOR."""

    def __init__(self, path: Path, batch_size: int = SQLITE_MULTISET_BATCH) -> None:
        if int(batch_size) <= 0:
            raise ValueError("SQLite multiset batch_size must be positive")
        self.db = sqlite3.connect(str(path))
        # 이 DB는 검증 도중에만 존재하고 종료 시 삭제된다. 중단되면 전체
        # trace를 처음부터 다시 계산하므로 rollback journal 내구성은 필요 없다.
        self.db.execute("PRAGMA journal_mode=OFF")
        self.db.execute("PRAGMA synchronous=OFF")
        self.db.execute("PRAGMA locking_mode=EXCLUSIVE")
        # worker별 bounded page cache; RAM에 DB 전체를 올리지는 않는다.
        self.db.execute(f"PRAGMA cache_size={-1024 * SQLITE_CACHE_MIB}")
        self.db.execute(
            "CREATE TABLE x(kind TEXT,digest BLOB,src INTEGER DEFAULT 0,out INTEGER DEFAULT 0,"
            "PRIMARY KEY(kind,digest)) WITHOUT ROWID"
        )
        self.total: dict[str, Counter[str]] = {"src": Counter(), "out": Counter()}
        self.batch_size = int(batch_size)
        self.pending: dict[str, list[tuple[str, bytes]]] = {"src": [], "out": []}

    def _flush(self, side: str) -> None:
        rows = self.pending[side]
        if not rows:
            return
        col = "src" if side == "src" else "out"
        sql = (
            f"INSERT INTO x(kind,digest,{col}) VALUES(?,?,1) "
            f"ON CONFLICT(kind,digest) DO UPDATE SET {col}={col}+1"
        )
        with self.db:
            self.db.executemany(sql, rows)
        rows.clear()

    def add(self, side: str, kind: str, digest: bytes) -> None:
        side = "src" if side == "src" else "out"
        self.pending[side].append((kind, sqlite3.Binary(digest)))
        self.total[side][kind] += 1
        if len(self.pending[side]) >= self.batch_size:
            self._flush(side)

    def result(self, max_examples: int) -> dict[str, Any]:
        self._flush("src")
        self._flush("out")
        mismatch: dict[str, Any] = {}
        for kind, digest, src, out in self.db.execute(
            "SELECT kind,digest,src,out FROM x WHERE src!=out ORDER BY kind,digest"
        ):
            item = mismatch.setdefault(
                kind, {"identity_count": 0, "source_only": 0, "zarr_only": 0, "examples": []}
            )
            item["identity_count"] += 1
            item["source_only"] += max(0, int(src) - int(out))
            item["zarr_only"] += max(0, int(out) - int(src))
            if len(item["examples"]) < max_examples:
                item["examples"].append(
                    {"sha256": bytes(digest).hex(), "source_count": int(src), "zarr_count": int(out)}
                )
        kinds = sorted(set(self.total["src"]) | set(self.total["out"]) | set(mismatch))
        return {
            "pass": not mismatch,
            "source_total": {key: int(self.total["src"][key]) for key in kinds},
            "zarr_total": {key: int(self.total["out"][key]) for key in kinds},
            "mismatches": mismatch,
        }

    def close(self) -> None:
        self.db.close()


def _make_bag(path,backend=None):
    backend=backend or os.environ.get('MFM_TRACE_MULTISET_BACKEND','auto')
    if backend=='auto':backend='partitioned'
    if backend=='sqlite':return _Bag(path)
    if backend!='partitioned':raise ValueError('MFM_TRACE_MULTISET_BACKEND must be partitioned|sqlite')
    return PartitionedBag(path,leaf_bytes=positive_env_int('MFM_TRACE_SORT_LEAF_MIB',64)<<20,
                          batch_records=SQLITE_MULTISET_BATCH)


def _cfg_bool(config: dict[str, Any], key: str, default: bool) -> bool:
    value = config.get(key, default)
    return value.strip().lower() in {"1", "true", "yes"} if isinstance(value, str) else bool(value)


def _reject(read: pysam.AlignedSegment, config: dict[str, Any]) -> str | None:
    checks = (
        ("unmapped", "exclude_unmapped", read.is_unmapped, True),
        ("secondary", "exclude_secondary", read.is_secondary, True),
        ("supplementary", "exclude_supplementary", read.is_supplementary, False),
        ("qcfail", "exclude_qcfail", read.is_qcfail, True),
        ("duplicate", "exclude_duplicate", read.is_duplicate, False),
        ("improper_pair", "require_proper_pair", read.is_paired and not read.is_proper_pair, False),
    )
    for reason, key, triggered, default in checks:
        if _cfg_bool(config, key, default) and triggered:
            return reason
    if int(read.mapping_quality) < int(config.get("min_mapping_quality", 30)):
        return "low_mapq"
    return "missing_sequence" if not read.query_sequence else None


def _mate(read: pysam.AlignedSegment) -> int:
    return 1 if read.is_read1 else (2 if read.is_read2 else 0)


def _source_key(read: pysam.AlignedSegment) -> tuple[Any, ...]:
    return (
        str(read.query_name), _mate(read), int(read.reference_start),
        int(read.reference_end or read.reference_start), read.cigarstring or "",
        bool(read.is_supplementary), bool(read.is_reverse),
    )


def _hp(read: pysam.AlignedSegment) -> int:
    try:
        value = int(read.get_tag("HP"))
    except (KeyError, TypeError, ValueError):
        return 0
    return value if value in (1, 2) else 0


def _ps(read: pysam.AlignedSegment) -> int:
    try:
        value = int(read.get_tag("PS"))
    except (KeyError, TypeError, ValueError):
        return 0
    return value if value > 0 else 0


def _sequence_and_cigar(read: pysam.AlignedSegment) -> tuple[str, int, tuple[int, ...]]:
    sequence = (read.query_sequence or "").upper()
    query_start = int(read.query_alignment_start or 0)
    query_end = int(read.query_alignment_end or len(sequence))
    packed = tuple(
        (int(length) << 4) | int(op)
        for op, length in (read.cigartuples or ())
        if int(op) not in (4, 5)
    )
    return sequence[query_start:query_end], query_start, packed


def _project(read: pysam.AlignedSegment, start: int, end: int) -> tuple[int, int, int, int] | None:
    """Independent identity-axis CIGAR projection; returns query/raw half-open ranges."""

    if not read.query_sequence or not read.cigartuples:
        return None
    q, r, previous = 0, int(read.reference_start), None
    qlo = qhi = rlo = rhi = None

    def choose(a: int, b: int, c: int, d: int) -> None:
        nonlocal qlo, qhi, rlo, rhi
        if b <= a:
            return
        qlo = a if qlo is None else min(qlo, a)
        qhi = b if qhi is None else max(qhi, b)
        rlo = c if rlo is None else min(rlo, c)
        rhi = d if rhi is None else max(rhi, d)

    for op, length0 in read.cigartuples:
        op, length = int(op), int(length0)
        if length <= 0:
            raise ValueError("invalid source CIGAR length")
        if op in (0, 7, 8):
            lo, hi = r, r + length
            overlap_lo, overlap_hi = max(lo, start), min(hi, end)
            if overlap_hi > overlap_lo:
                choose(q + overlap_lo - lo, q + overlap_hi - lo, overlap_lo, overlap_hi)
            q, r, previous = q + length, hi, hi - 1
        elif op == 1:
            anchor = previous + 1 if previous is not None else r
            if start <= anchor < end:
                choose(q, q + length, anchor, anchor + 1)
            q += length
        elif op in (2, 3):
            r += length
        elif op == 4:
            q += length
        elif op in (5, 6):
            pass
        else:
            raise ValueError(f"unsupported source CIGAR op: {op}")
    if None in (qlo, qhi, rlo, rhi):
        return None
    return int(qlo), int(qhi), int(rlo), int(rhi)


def _match_intervals(ref_start: int, cigar: Iterable[int]) -> list[tuple[int, int]]:
    pos, intervals = int(ref_start), []
    for value in cigar:
        op, length = int(value) & 15, int(value) >> 4
        if op in (0, 7, 8):
            intervals.append((pos, pos + length))
            pos += length
        elif op in (2, 3):
            pos += length
        elif op in (1, 6):
            continue
        else:
            raise ValueError(f"unsupported packed CIGAR op: {op}")
    return intervals


def _merge(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    for start, end in sorted((int(a), int(b)) for a, b in intervals if b > a):
        if not result or start > result[-1][1]:
            result.append((start, end))
        else:
            result[-1] = (result[-1][0], max(result[-1][1], end))
    return result


def _stage05_windows(path: Path, contig: str, hap: str) -> list[tuple[int, int]]:
    rows: list[tuple[int, int, bool]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"contig", "haplotype", "raw_start", "actual_length", "kept"}
        if required - set(reader.fieldnames or ()):
            raise ValueError("Stage 05 window index schema mismatch")
        for row in reader:
            if row.get("contig") == contig and row.get("haplotype") == hap:
                kept = str(row["kept"]).strip().lower()
                if kept not in {"0", "1", "false", "true"}:
                    raise ValueError("invalid Stage 05 kept value")
                rows.append((int(row["raw_start"]), int(row["actual_length"]), kept in {"1", "true"}))
    keys = [(a, b) for a, b, _ in rows]
    if not rows or len(keys) != len(set(keys)):
        raise ValueError(f"invalid/duplicate Stage 05 windows for {hap}/{contig}")
    kept = [(a, b) for a, b, c in rows if c]
    if not kept or kept != sorted(kept):
        raise ValueError(f"empty/unsorted kept Stage 05 windows for {hap}/{contig}")
    return kept


def _windows(hap_group: zarr.Group) -> list[tuple[int, int, int]]:
    starts = np.asarray(hap_group["raw_start"][:], dtype=np.int64)
    actual = np.asarray(hap_group["actual_length"][:], dtype=np.int64)
    stored = np.asarray(hap_group["lengths"][:], dtype=np.int64)
    if starts.shape != actual.shape or starts.shape != stored.shape:
        raise ValueError("HAP window array shape mismatch")
    return [(int(a), int(b), int(c)) for a, b, c in zip(starts, actual, stored, strict=True)]


def _check_hap_windows(
    group: zarr.Group, fasta_path: Path, stage05: list[tuple[int, int]], contig: str, max_examples: int
) -> tuple[list[tuple[int, int, int]], dict[str, Any]]:
    windows = _windows(group)
    observed_keys = [(a, b) for a, b, _ in windows]
    selection_ok = observed_keys == stage05
    offsets = np.asarray(group["offsets"][:], dtype=np.int64)
    if offsets.size != len(windows) + 1:
        raise ValueError("HAP packed offsets mismatch")
    bad, examples = 0, []
    with pysam.FastaFile(str(fasta_path)) as fasta:
        if contig not in fasta.references:
            raise ValueError("HAP FASTA contig missing")
        for index, (start, actual, stored) in enumerate(windows):
            observed = _unpack(group["bases_packed"][offsets[index]:offsets[index + 1]], stored, HAP_BASES, "hap")
            expected = fasta.fetch(contig, start, start + actual).upper() + "P" * (stored - actual)
            if observed != expected:
                bad += 1
                if len(examples) < max_examples:
                    examples.append({"window_index": index, "raw_start": start})
    return windows, {
        "pass": selection_ok and bad == 0,
        "stage05_selection_pass": selection_ok,
        "hap_fasta_mismatch_windows": bad,
        "examples": examples,
    }


def _source(
    alignment: Path, target: Path, contig: str, windows: list[tuple[int, int, int]], config: dict[str, Any],
    records: ReplayRecords | None = None,
) -> tuple[Iterable[dict[str, Any]], dict[str, int]]:
    audit, seen = Counter(), {}
    intervals = _merge((start, start + actual) for start, actual, _ in windows)
    lookup = WindowLookup(windows)
    heartbeat = time.monotonic()
    pending=[];pending_bytes=0
    def flush_hashes():
        if not pending:return
        hashed=source_digests([item[0] for item in pending])
        for (_canonical,record),digest in zip(pending,hashed):
            record['digest']=bytes(digest)
            if records is not None:records.append(record)
        pending.clear()
    with pysam.AlignmentFile(str(alignment), "rb", reference_filename=str(target)) as bam:
        name = contig if contig in bam.references else contig.removeprefix("chr")
        if name not in bam.references:
            raise ValueError("source CRAM contig missing")
        for start, end in intervals:
            for read in bam.fetch(name, start, end):
                audit["fetched_records"] += 1
                if audit["fetched_records"] % 100000 == 0 and time.monotonic() - heartbeat >= 60:
                    print(f"[trace] source contig={contig} alignment={alignment.name} "
                          f"fetched={audit['fetched_records']} retained={audit['retained_alignments']}", flush=True)
                    heartbeat = time.monotonic()
                if records is not None and audit["fetched_records"] % 10000 == 0:
                    # Indexed coordinate fetch에서 이전 위치로 다시 들어오는 정렬은
                    # 다음 interval에도 닿는 read뿐이다. 경계의 trailing insertion은 유지한다.
                    boundary = int(read.reference_start)
                    expired = [key for key, end in seen.items() if end < boundary]
                    for key in expired:
                        del seen[key]
                reason = _reject(read, config)
                if reason:
                    audit[f"excluded_{reason}"] += 1
                    continue
                key = _source_key(read)
                if key in seen:
                    audit["deduplicated_repeated_alignment"] += 1
                    continue
                sequence, clip_start, cigar = _sequence_and_cigar(read)
                if not sequence:
                    audit["excluded_empty_mapped_query"] += 1
                    continue
                segments = []
                for index in lookup.candidates(read):
                    wstart, actual, _stored = windows[index]
                    project = _project(read, wstart, wstart + actual)
                    if project is None:
                        continue
                    q0, q1, r0, r1 = project
                    q0, q1, h0, h1 = q0 - clip_start, q1 - clip_start, r0 - wstart, r1 - wstart
                    if q0 < 0 or q1 <= q0 or q1 > len(sequence) or h0 < 0 or h1 <= h0 or h1 > actual:
                        segments = []
                        audit["excluded_invalid_segment_coordinates"] += 1
                        break
                    segments.append((index, q0, q1, h0, h1))
                if not segments:
                    audit["excluded_unprojectable_or_outside_windows"] += 1
                    continue
                canonical = (
                    sequence, _hp(read), _ps(read), int(read.reference_start),
                    int(read.reference_end or read.reference_start), int(read.is_supplementary), cigar,
                )
                record = {
                    "segment_type": int(read.is_supplementary),
                    "digest": None, "segments": segments, "cigar": cigar,
                    "ref_start": int(read.reference_start), "strand": -1 if read.is_reverse else 1,
                }
                if records is None:
                    seen[key] = record
                else:
                    seen[key] = int(read.reference_end or read.reference_start)
                pending.append((canonical,record));pending_bytes+=len(sequence)+8*len(cigar)
                if len(pending)>=4096 or pending_bytes>=16<<20:
                    flush_hashes();pending_bytes=0
                audit["retained_alignments"] += 1
                if read.is_supplementary:
                    audit["retained_supplementary"] += 1
                if read.is_duplicate:
                    audit["retained_duplicate"] += 1
    flush_hashes()
    return (list(seen.values()) if records is None else records), {key: int(value) for key, value in sorted(audit.items())}


def _output_reads(root: zarr.Group, bag: _Bag, prefetch_records: int = ZARR_PREFETCH_RECORDS) -> np.ndarray:
    """Zarr metadata/CIGAR를 연속 범위로 읽고 기존 canonical digest 순서를 보존한다."""
    if int(prefetch_records) <= 0:
        raise ValueError("prefetch_records must be positive")
    reads, cigar = root["reads"], root["cigar"]["ops_packed"]
    count = int(reads["lengths"].shape[0])
    if int(reads["offsets"].shape[0]) != count + 1:
        raise ValueError("output read offsets mismatch")
    result = np.empty(count, dtype="V32")
    decoded_cache=DecodedChunkCache(nonnegative_env_int("MFM_TRACE_DECODED_CACHE_MIB",64)<<20)
    metadata_names = ("lengths", "hp", "ps", "ref_start", "ref_end", "segment_type", "cigar_offset", "cigar_len")
    max_packed=positive_env_int("MFM_TRACE_PREFETCH_PACKED_MIB",16)<<20
    begin=0
    while begin<count:
        end = min(count, begin + int(prefetch_records))
        offsets = np.asarray(reads["offsets"][begin:end + 1], dtype=np.int64)
        if offsets.size!=end-begin+1 or np.any(offsets[1:]<offsets[:-1]):raise ValueError('output read offsets mismatch')
        if int(offsets[-1])-int(offsets[0])>max_packed:
            if int(offsets[1])-int(offsets[0])>max_packed:raise ValueError('trace single read exceeds packed batch budget')
            selected=max(1,int(np.searchsorted(offsets,int(offsets[0])+max_packed,side='right'))-1)
            end=begin+selected;offsets=offsets[:selected+1]
        metadata = {name: np.asarray(reads[name][begin:end]) for name in metadata_names}
        if offsets.size != end - begin + 1 or any(values.size != end - begin for values in metadata.values()):
            raise ValueError("output read metadata batch shape mismatch")
        packed_begin, packed_end = int(offsets[0]), int(offsets[-1])
        packed_array=reads["bases_packed"]
        packed=np.asarray(decoded_cache.read(packed_array,packed_begin,packed_end) if isinstance(packed_array,zarr.Array) else packed_array[packed_begin:packed_end],dtype=np.uint8)
        cigar_offsets = np.asarray(metadata["cigar_offset"], dtype=np.int64)
        cigar_lengths = np.asarray(metadata["cigar_len"], dtype=np.int64)
        if cigar_lengths.size and np.any(cigar_lengths < 0):
            raise ValueError("negative output CIGAR length")
        cigar_begin = int(cigar_offsets.min()) if cigar_offsets.size else 0
        cigar_end = int(np.max(cigar_offsets + cigar_lengths)) if cigar_offsets.size else 0
        cigar_values=np.asarray(decoded_cache.read(cigar,cigar_begin,cigar_end) if isinstance(cigar,zarr.Array) else cigar[cigar_begin:cigar_end])
        columns=np.column_stack((metadata['hp'],metadata['ps'],metadata['ref_start'],metadata['ref_end'],
            metadata['segment_type'],cigar_offsets-cigar_begin,cigar_lengths,metadata['lengths']))
        batch=digests(packed,offsets-packed_begin,columns,cigar_values,packed=True)
        result[begin:end]=batch
        if hasattr(bag,'add_many'):bag.add_many('out','read',batch)
        else:
            for digest in batch:bag.add('out','read',bytes(digest))
        begin=end
    print('[TRACE-CACHE] ' + json.dumps({
        'limit_bytes': decoded_cache.limit, 'retained_bytes': decoded_cache.bytes,
        'hits': decoded_cache.hits, 'misses': decoded_cache.misses,
        'packed_batch_bytes': max_packed, 'records': count,
    }), flush=True)
    return result


def _seg_digest(read: bytes, index: int, q0: int, q1: int, h0: int, h1: int, strand: int) -> bytes:
    return segment_digest(read,index,q0,q1,h0,h1,strand)


def _check_projection(
    segments: zarr.Group, depth_group: zarr.Group, windows: list[tuple[int, int, int]],
    source: Iterable[dict[str, Any]], read_ids: Iterable[bytes], bag: _Bag, max_examples: int, allowed_segment_types: tuple[int, ...],
    *, check_segments: bool = True,
) -> dict[str, Any]:
    expected_depth = DepthAccumulator(windows)
    source_hashes=SegmentBatch(bag)
    allowed = {int(value) for value in allowed_segment_types}
    for record in source:
        include_depth = int(record["segment_type"]) in allowed
        intervals = _match_intervals(record["ref_start"], record["cigar"]) if include_depth or check_segments else []
        for index, q0, q1, h0, h1 in record["segments"]:
            if check_segments:
                source_hashes.add(record["digest"], index, q0, q1, h0, h1, record["strand"])
            wstart, actual, _stored = windows[index]
            if not include_depth:
                continue
            for r0, r1 in intervals:
                a, b = max(r0, wstart), min(r1, wstart + actual)
                if b > a:
                    expected_depth.add(index, a - wstart, b - wstart)
    source_hashes.flush()
    soff = np.asarray(segments["window_read_offsets"][:], dtype=np.int64)
    doff = np.asarray(depth_group["offsets"][:], dtype=np.int64)
    if soff.size != len(windows) + 1 or doff.size != len(windows) + 1:
        raise ValueError("output segment/depth offsets mismatch")
    bad_windows = bad_bases = 0
    examples = []
    segment_names = ("read_id", "query_start", "query_end", "hap_start", "hap_end", "strand")
    for index, (wstart, _actual, _stored) in enumerate(windows):
        segment_begin, segment_end = int(soff[index]), int(soff[index + 1])
        segment_values = {
            name: np.asarray(segments[name][segment_begin:segment_end])
            for name in segment_names
        } if check_segments else {}
        segment_count = segment_end - segment_begin
        if any(values.size != segment_count for values in segment_values.values()):
            raise ValueError("output segment metadata batch shape mismatch")
        for begin in range(0,segment_count if check_segments else 0,8192):
            stop=min(segment_count,begin+8192)
            ids=np.asarray(segment_values['read_id'][begin:stop],dtype=np.int64)
            if np.any(ids<0) or np.any(ids>=len(read_ids)):raise ValueError('invalid output segment read_id')
            fields=np.column_stack((np.full(len(ids),index,dtype=np.int64),
                *(segment_values[name][begin:stop] for name in segment_names[1:])))
            hashed=segment_digests(np.asarray(read_ids,dtype='V32')[ids],fields)
            add_hashes(bag,'out','segment',hashed)
        observed = np.asarray(depth_group["depth"][doff[index]:doff[index + 1]], dtype=np.uint32)
        expected = expected_depth.finish_window(index)
        if observed.shape != expected.shape:
            raise ValueError("output depth length mismatch")
        differing = np.flatnonzero(observed != expected)
        if differing.size:
            bad_windows += 1
            bad_bases += int(differing.size)
            if len(examples) < max_examples:
                base = int(differing[0])
                examples.append({
                    "window_index": index, "raw_position": wstart + base,
                    "source_depth": int(expected[base]), "zarr_depth": int(observed[base]),
                })
    return {"depth_mismatch_windows": bad_windows, "depth_mismatch_bases": bad_bases, "examples": examples}


def _paths(metadata: dict[str, Any], hap: str) -> tuple[Path, Path, Path, Path]:
    provenance = metadata.get("input_provenance", {})
    data = provenance.get("haplotypes", {}).get(hap, {})
    try:
        alignment = Path(data["alignment"]["path"])
        target = Path(data["target_fasta"]["path"])
        fasta = Path(data["hap_fasta"]["path"])
        index = Path(provenance["window_index"]["path"])
    except KeyError as exc:
        raise ValueError(f"malformed provenance for {hap}") from exc
    # canonical CRAM과 동일한 원래 링크를 열어 링크 옆 .crai를 찾음.
    logical=metadata.get('source_alignment',{})
    if isinstance(logical,dict) and logical.get(hap):
        candidate=Path(logical[hap])
        if candidate.exists():
            if not candidate.samefile(alignment):raise ValueError(f'trace alignment provenance mismatch for {hap}')
            alignment=candidate
    for path in (alignment, target, fasta, index):
        if not path.exists():
            raise FileNotFoundError(f"trace input missing: {path}")
    return alignment, target, fasta, index


def _phase_done(contig: str, hap: str, phase: str, started: float) -> float:
    now = time.monotonic()
    print(f"[trace] phase contig={contig} hap={hap} name={phase} elapsed_s={now - started:.3f}", flush=True)
    return now


def _validate_hap(path: Path, metadata: dict[str, Any], hap: str, scratch: Path, max_examples: int) -> dict[str, Any]:
    phase_started = time.monotonic()
    alignment, target, fasta, index = _paths(metadata, hap)
    config = metadata["input_provenance"].get("export_config")
    if not isinstance(config, dict):
        raise ValueError("missing provenance export_config")
    primary_segment_types = tuple(int(value) for value in config.get("depth_segment_types", (0, 1)))
    supplementary_segment_types = tuple(int(value) for value in config.get("supplementary_depth_segment_types", ()))
    contig = str(metadata["contig"])
    hap_root = zarr.open_group(str(path / "hap_sequence.zarr"), mode="r")
    readset = zarr.open_group(str(path / f"readset.{hap}.zarr"), mode="r")
    depth_root = zarr.open_group(str(path / "depth.zarr"), mode="r")
    stage05 = _stage05_windows(index, contig, hap)
    max_windows = metadata["input_provenance"].get("max_windows")
    if max_windows is not None:
        if isinstance(max_windows, bool) or int(max_windows) <= 0:
            raise ValueError("invalid provenance max_windows")
        stage05 = stage05[:int(max_windows)]
    windows, window_qc = _check_hap_windows(hap_root[hap], fasta, stage05, contig, max_examples)
    phase_started = _phase_done(contig, hap, "hap_windows", phase_started)
    with tempfile.TemporaryDirectory(prefix=f".stage70_trace_{hap}_", dir=str(scratch)) as tmp:
        records = ReplayRecords(Path(tmp) / "source.records")
        bag = _make_bag(Path(tmp) / "multiset.sqlite3",metadata.get("_trace_multiset_backend"))
        try:
            source, audit = _source(alignment, target, contig, windows, config, records)
            phase_started = _phase_done(contig, hap, "source_projection", phase_started)
            for record in source:
                bag.add("src", "read", record["digest"])
            read_ids = _output_reads(readset, bag)
            phase_started = _phase_done(contig, hap, "read_multiset", phase_started)
            projection = _check_projection(readset["segments"][hap], depth_root[hap], windows, source, read_ids, bag, max_examples, primary_segment_types)
            phase_started = _phase_done(contig, hap, "segment_primary_depth", phase_started)
            multiset = bag.result(max_examples)
            phase_started = _phase_done(contig, hap, "multiset_compare", phase_started)
            supplementary_projection = None
            if supplementary_segment_types:
                supplementary_path = path / "supplementary_depth.zarr"
                if not supplementary_path.exists():
                    raise ValueError("supplementary_depth.zarr is required by export_config")
                supplementary_root = zarr.open_group(str(supplementary_path), mode="r")
                # 첫 pass에서 모든 primary/supplementary segment를 SHA multiset으로
                # 대조했다. 두 번째는 supplementary depth만 재구성한다.
                supplementary_projection = _check_projection(
                    readset["segments"][hap], supplementary_root[hap], windows, source,
                    read_ids, bag, max_examples, supplementary_segment_types, check_segments=False,
                )
                _phase_done(contig, hap, "supplementary_depth", phase_started)
        finally:
            records.close()
            bag.close()
    passed = window_qc["pass"] and multiset["pass"] and projection["depth_mismatch_bases"] == 0 and (supplementary_projection is None or supplementary_projection["depth_mismatch_bases"] == 0)
    return {
        "multiset_execution": getattr(bag,"stats",{"backend":"sqlite"}),
        "canonical_hash_backend":hash_backend(),
        "supplementary_projection": supplementary_projection,
        "status": "pass" if passed else "fail",
        "source": {"alignment": str(alignment), "target_fasta": str(target), "hap_fasta": str(fasta), "filter_config": config},
        "windows": window_qc, "filter_audit": audit, "record_multiset": multiset, "projection": projection,
    }


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def validate_stage70_trace_contig(
    compact_contig: Path, report_path: Path | None = None, scratch_dir: Path | None = None,
    max_examples: int = 20, workers: int = 1,
) -> dict[str, Any]:
    """Read-only source-to-Zarr validation with optional HAP-level workers."""

    if int(workers) <= 0:
        raise ValueError("trace workers must be positive")
    compact_contig = Path(compact_contig)
    report, metadata = _prepare_contig_report(compact_contig)
    if metadata is not None:
        scratch = Path(scratch_dir) if scratch_dir else compact_contig.parent.parent / ".trace_scratch"
        scratch.mkdir(parents=True, exist_ok=True)
        key = str(metadata.get("contig", compact_contig.name))
        tasks = [
            ((key, hap), compact_contig, metadata, hap, scratch, int(max_examples))
            for hap in HAPS
        ]
        hap_results = _run_hap_tasks(tasks, min(int(workers), len(HAPS)))
        report["by_haplotype"] = {hap: hap_results[(key, hap)] for hap in HAPS}
        report["status"] = (
            "pass" if all(item.get("status") == "pass" for item in report["by_haplotype"].values())
            else "fail"
        )
    if report_path:
        _write_json(Path(report_path), report)
    return report


def _prepare_contig_report(compact_contig: Path) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """Load the immutable task context; malformed contigs remain isolated failures."""

    report: dict[str, Any] = {
        "schema": TRACE_SCHEMA, "compact_contig": str(compact_contig), "status": "fail",
        "by_haplotype": {}, "errors": [],
    }
    try:
        metadata = json.loads((compact_contig / "contig_metadata.json").read_text(encoding="utf-8"))
        if metadata.get("read_coordinate_mode") != "personal_realigned" or metadata.get("readset_layout") != "hap_split":
            raise ValueError("only personal_realigned hap_split compact is supported")
        if not metadata.get("input_provenance", {}).get("fingerprint_sha256"):
            raise ValueError("compact input provenance is required")
        report.update({
            "sample_id": str(metadata.get("sample_id", "")), "contig": str(metadata.get("contig", "")),
            "input_provenance_fingerprint": metadata["input_provenance"]["fingerprint_sha256"],
        })
        return report, metadata
    except Exception as exc:
        report["errors"].append(f"{type(exc).__name__}: {exc}")
        return report, None


def _safe_validate_hap(
    compact_contig: Path, metadata: dict[str, Any], hap: str, scratch: Path, max_examples: int,
) -> dict[str, Any]:
    """Process boundary: convert one HAP failure to structured QC."""

    started = time.monotonic()
    label = f"contig={metadata.get('contig', compact_contig.name)} hap={hap} pid={os.getpid()}"
    print(f"[trace] task start {label}", flush=True)
    try:
        from .pipeline.checkpoints import run_step, file_identity, tree_identity, digest
        checkpoint_dir = metadata.get("_trace_checkpoint_dir")
        if checkpoint_dir:
            inputs = list(_paths(metadata, hap))
            for base in list(inputs):
                for suffix in (".fai", ".crai", ".bai"):
                    extra = Path(str(base) + suffix)
                    if extra.is_file(): inputs.append(extra)
            module = Path(__file__).parent
            code = {name: digest((module / name).read_text()) for name in
                    ("stage70_trace_validator.py", "trace_acceleration.py", "compact/trace_store.py", "compact/trace_decode.py", "compact/trace_multiset.py", "compact/trace_resources.py", "window_chunk_cache.py", "pipeline/checkpoints.py")}
            def identity():
                return {"inputs": [file_identity(p) for p in inputs], "compact": tree_identity(compact_contig)}
            result = run_step(Path(checkpoint_dir) / str(metadata["contig"]) / (hap + ".json"),
                              {"metadata": metadata, "max_examples": max_examples, "code": code},
                              lambda: _validate_hap(compact_contig, metadata, hap, scratch, max_examples),
                              lambda: tree_identity(compact_contig), inputs_identity=identity)
        else:
            result = _validate_hap(compact_contig, metadata, hap, scratch, max_examples)
    except Exception as exc:
        result = {"status": "fail", "error": f"{type(exc).__name__}: {exc}"}
    print(f"[trace] task done {label} status={result['status']} elapsed_s={time.monotonic() - started:.3f}", flush=True)
    return result


def _run_hap_tasks(
    tasks: list[tuple[tuple[str, str], Path, dict[str, Any], str, Path, int]], workers: int,
) -> dict[tuple[str, str], dict[str, Any]]:
    """Run one flat (contig,HAP) pool so nested pools cannot oversubscribe."""

    if not tasks:
        return {}
    estimates = []
    for key, path, metadata, hap, scratch, _examples in tasks:
        try:
            group = zarr.open_group(str(path / f"readset.{hap}.zarr"), mode="r")
            hap_group = zarr.open_group(str(path / "hap_sequence.zarr"), mode="r")[hap]
            n_reads = int(group["reads"]["lengths"].shape[0])
            bases = int(hap_group["offsets"][-1]) * 2
            segments=int(group["segments"][hap]["read_id"].shape[0]);cigar_ops=int(group["cigar/ops_packed"].shape[0])
            sort_extra=(8*positive_env_int('MFM_TRACE_SORT_LEAF_MIB',64)+nonnegative_env_int('MFM_TRACE_DECODED_CACHE_MIB',64)+16*positive_env_int('MFM_TRACE_PREFETCH_PACKED_MIB',16))*2**20
            peak = int((2 * GIB + n_reads * 32 + bases * 4+sort_extra) * 1.3)
            ceiling = max(8 * GIB, ((peak + 2 * GIB + (256 << 20) - 1) // (256 << 20)) * (256 << 20))
        except Exception:
            # 개별 입력 오류는 기존 task의 구조화된 QC 실패 경로에서 처리한다.
            n_reads = 0
            segments=0;cigar_ops=0
            peak = 6 * GIB
            ceiling = 8 * GIB
        estimates.append({"contig": key[0], "hap": hap, "reads": n_reads,
                          "estimated_peak_bytes": peak, "worker_ceiling_bytes": ceiling,**scratch_estimate(n_reads,segments,cigar_ops)})
    elastic = elastic_cpu_enabled("MFM_TRACE_JOBS_REQUESTED",thread_cost_accounted=True)
    scratch_plans={}
    for scratch in set(task[4] for task in tasks):
        chosen=[e for e,t in zip(estimates,tasks) if t[4]==scratch]
        plan=select_backend(scratch,chosen,os.environ.get('MFM_TRACE_MULTISET_BACKEND','auto'),workers=workers,elastic=elastic);scratch_plans[scratch]=plan
        print(f"[trace] scratch backend={plan['backend']} concurrent_cap={plan['max_concurrent_workers']} peak_bytes={plan['estimated_scratch_bytes']} all_tasks_bytes={plan['all_tasks_scratch_bytes']} free_bytes={plan['free_bytes']}",flush=True)
        if plan['memory_filesystem']:
            for e in chosen:e['worker_ceiling_bytes']+=e[plan['backend']+'_estimate_bytes']
    tasks=[(key,path,{**meta,'_trace_multiset_backend':scratch_plans[scratch]['backend']},hap,scratch,examples)
           for key,path,meta,hap,scratch,examples in tasks]
    ceiling = max(row["worker_ceiling_bytes"] for row in estimates)
    if ceiling + 4 * GIB > memory_limit() * 9 // 10:
        raise ValueError("독립 검증 한 작업의 RAM 예측량이 전체 예산을 초과합니다.")
    for scratch in set(task[4] for task in tasks):
        _write_json(scratch / "memory_preflight_trace.json", {
            "schema": "mfm_trace_memory_preflight_v1", "memory_limit_bytes": memory_limit(),
            "worker_ceiling_bytes": ceiling, "tasks": estimates,"scratch":scratch_plans[scratch],
            "source_records": "디스크 재생 저장; 출력 digest는 read당 32바이트 고정 폭"})
    results: dict[tuple[str, str], dict[str, Any]] = {}
    for index, result, error in cpu_results(
        _safe_validate_hap, [tuple(task[1:]) for task in tasks], workers,
        elastic=elastic,
        max_workers=min(plan["max_concurrent_workers"] for plan in scratch_plans.values()),
        mp_context=get_context("spawn"),
        memory_per_worker_bytes=ceiling,
        cpus_per_worker=TRACE_BLOSC_THREADS,
    ):
        key = tasks[index][0]
        results[key] = ({"status": "fail", "error": f"{type(error).__name__}: {error}"}
                        if error is not None else result)
        print(f"[trace] progress completed={len(results)}/{len(tasks)} contig={key[0]} hap={key[1]} status={results[key]['status']}", flush=True)
    return results


def validate_stage70_trace_root(
    compact_root: Path, contigs: Iterable[str], report_root: Path, max_examples: int = 20,
    scratch_dir: Path | None = None, workers: int = 1,
) -> dict[str, Any]:
    """Validate a flat contig×HAP task pool and materialize reports in input order."""

    if int(workers) <= 0:
        raise ValueError("trace workers must be positive")
    compact_root, report_root = Path(compact_root), Path(report_root)
    if compact_root.resolve() == report_root.resolve() or compact_root.resolve() in report_root.resolve().parents:
        raise ValueError("report_root must be outside compact_root")
    report_root.mkdir(parents=True, exist_ok=True)
    scratch_root = Path(scratch_dir) if scratch_dir is not None else report_root
    scratch_root.mkdir(parents=True, exist_ok=True)

    # The old dict result already collapsed duplicate contigs. Make that behavior
    # explicit while retaining the caller's first-seen order.
    contig_order = list(dict.fromkeys(str(contig) for contig in contigs))
    contexts: dict[str, tuple[dict[str, Any], dict[str, Any] | None]] = {}
    tasks: list[tuple[tuple[str, str], Path, dict[str, Any], str, Path, int]] = []
    for contig in contig_order:
        compact_contig = compact_root / contig
        report, metadata = _prepare_contig_report(compact_contig)
        if metadata is not None:
            metadata = {**metadata, "_trace_checkpoint_dir": str(report_root / "hap_steps")}
        contexts[contig] = (report, metadata)
        if metadata is not None:
            tasks.extend(
                ((contig, hap), compact_contig, metadata, hap, scratch_root, int(max_examples))
                for hap in HAPS
            )

    task_results = _run_hap_tasks(tasks, int(workers))
    results: dict[str, dict[str, Any]] = {}
    for contig in contig_order:
        report, metadata = contexts[contig]
        if metadata is not None:
            report["by_haplotype"] = {hap: task_results[(contig, hap)] for hap in HAPS}
            report["status"] = (
                "pass" if all(item.get("status") == "pass" for item in report["by_haplotype"].values())
                else "fail"
            )
        results[contig] = report
        _write_json(report_root / f"{contig}.trace_qc.json", report)

    summary = {
        "schema": TRACE_SCHEMA, "compact_root": str(compact_root),
        "status": "pass" if results and all(item["status"] == "pass" for item in results.values()) else "fail",
        "contigs": {contig: results[contig]["status"] for contig in contig_order},
        "resources": {
            "requested_workers": int(workers),
            "effective_workers": min(int(workers), len(tasks)),
            "task_count": len(tasks),
            "blosc_threads_per_worker": TRACE_BLOSC_THREADS,
        },
    }
    _write_json(report_root / "trace_manifest.json", summary)
    return {"summary": summary, "results": results}
