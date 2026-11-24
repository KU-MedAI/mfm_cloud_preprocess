
from __future__ import annotations

import json
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import zarr
from numcodecs import Blosc

from .timing import StepTimer

SEQUENCE_ALLOWED = frozenset({"A", "C", "G", "T", "N", "P"})
READ_ALLOWED = frozenset({"A", "C", "G", "T", "N"})
SYMBOL_ORDER = ("A", "C", "G", "T", "N", "P")
HAP_ATTR_COPY_KEYS = (
    "coordinate_system",
    "coordinate_map_source",
    "source_chain_path",
    "source_individual_fasta",
    "source_bam_path",
    "random_offset",
    "random_offset_unit",
    "valid_window_start_index",
    "sequence_valid_window_count_total",
    "coordinate_eligible_window_count_total",
    "coordinate_excluded_window_count",
    "coordinate_excluded_bases",
    "coordinate_coverage_threshold",
)

ASCII_TO_CODE = np.full(256, 255, dtype=np.uint8)
for base, code in {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4, "P": 5}.items():
    ASCII_TO_CODE[ord(base)] = code
    ASCII_TO_CODE[ord(base.lower())] = code


@dataclass
class LegacyConvertConfig:
    sample_id: str
    source_sample_id: str
    legacy_root: Path
    output_sample_root: Path
    dataset_role: str
    read_type: str
    project_id: str
    contigs: list[str] | None = None
    overwrite: bool = False
    chunk_records: int = 1024
    hap_invalid_policy: str = "to_n"
    rna_bigwigs: list[Path] = field(default_factory=list)
    source_wgs: list[Path] = field(default_factory=list)
    source_vcfs: list[Path] = field(default_factory=list)
    layout_name: str = "legacy_shared_hap_4bit_acgtnp"


def run_legacy_convert(config: LegacyConvertConfig) -> dict[str, Any]:
    timer = StepTimer()
    with timer.step("prepare_output_root"):
        prepare_sample_root(config)

    contigs = resolve_contigs(config.legacy_root, config.contigs)
    contig_summaries = []
    for contig in contigs:
        with timer.step(f"convert_{contig}"):
            contig_summaries.append(convert_contig(config, contig))

    with timer.step("write_sample_metadata"):
        summary = {
            "sample_id": config.sample_id,
            "source_sample_id": config.source_sample_id,
            "dataset_role": config.dataset_role,
            "read_type": config.read_type,
            "project_id": config.project_id,
            "legacy_root": str(config.legacy_root),
            "output_sample_root": str(config.output_sample_root),
            "encoding": "4bit_acgtnp",
            "hap_invalid_policy": config.hap_invalid_policy,
            "sequence_allowed_symbols": format_symbols(SEQUENCE_ALLOWED),
            "read_allowed_symbols": format_symbols(READ_ALLOWED),
            "hap_invalid_replacement": "N" if config.hap_invalid_policy == "to_n" else "",
            "contigs": contig_summaries,
            "timing": timer.records,
        }
        write_metadata(config, summary)
    return summary


def prepare_sample_root(config: LegacyConvertConfig) -> None:
    config.output_sample_root.mkdir(parents=True, exist_ok=True)
    for name in ("metadata", "compact", "tracks", "qc", "logs"):
        (config.output_sample_root / name).mkdir(parents=True, exist_ok=True)


def resolve_contigs(legacy_root: Path, requested: list[str] | None) -> list[str]:
    available = sorted(path.stem for path in legacy_root.glob("chr*.zarr"))
    if requested:
        missing = [contig for contig in requested if not (legacy_root / f"{contig}.zarr").exists()]
        if missing:
            raise FileNotFoundError(f"requested contigs missing in {legacy_root}: {','.join(missing)}")
        return requested
    return available


def convert_contig(config: LegacyConvertConfig, contig: str) -> dict[str, Any]:
    source = config.legacy_root / f"{contig}.zarr"
    target_root = config.output_sample_root / "compact" / config.layout_name / contig
    if target_root.exists():
        if not config.overwrite:
            return {
                "contig": contig,
                "status": "skipped_existing",
                "source": str(source),
                "target": str(target_root),
            }
        shutil.rmtree(target_root)
    target_root.mkdir(parents=True, exist_ok=True)

    src = zarr.open_group(str(source), mode="r")
    hap_seq_root = zarr.open_group(str(target_root / "hap_sequence.zarr"), mode="w")
    readset_root = zarr.open_group(str(target_root / "readset.zarr"), mode="w")
    comp = Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE)

    hap_windows = []
    hap_summaries: dict[str, Any] = {}
    for hap in ("HAP1", "HAP2"):
        hap_summary = convert_hap_group(
            src[hap],
            hap_seq_root.create_group(hap),
            readset_root,
            hap,
            comp,
            config.chunk_records,
            config.hap_invalid_policy,
        )
        hap_summaries[hap] = hap_summary
        hap_windows.extend(make_window_rows(contig, hap, hap_summary))

    read_summary = convert_reads_group(src, src["reads"], readset_root, readset_root.create_group("reads"), comp, config.chunk_records)
    write_windows_tsv(target_root / "windows.tsv", hap_windows)
    write_contig_metadata(target_root, config, contig, source, hap_windows, hap_summaries, read_summary)

    size_bytes = directory_size(target_root)
    return {
        "contig": contig,
        "status": "converted",
        "source": str(source),
        "target": str(target_root),
        "window_rows": len(hap_windows),
        "read_count": read_summary["read_count"],
        "read_bases": read_summary["total_bases"],
        "hap_invalid_base_count": sum(int(item.get("invalid_base_count", 0)) for item in hap_summaries.values()),
        "hap_invalid_base_symbols": merge_count_dicts([item.get("invalid_base_symbols", {}) for item in hap_summaries.values()]),
        "size_bytes": size_bytes,
        "size_human": human_bytes(size_bytes),
    }


def convert_hap_group(src_hap: zarr.Group, dst_hap: zarr.Group, readset_root: zarr.Group, hap: str, comp: Blosc, chunk_records: int, invalid_policy: str) -> dict[str, Any]:
    copied_attrs = selected_hap_attrs(src_hap)
    dst_hap.attrs.update({
        "encoding": "4bit_acgtnp",
        "symbols": format_symbols(SEQUENCE_ALLOWED),
        "kind": "hap_window_sequence",
        "invalid_base_policy": invalid_policy,
        "invalid_base_replacement": "N" if invalid_policy == "to_n" else "",
        **copied_attrs,
    })
    seq_summary = pack_object_sequence_array(src_hap["sequence"], dst_hap, SEQUENCE_ALLOWED, comp, chunk_records, invalid_policy)
    copy_array(src_hap, dst_hap, "raw_start", "raw_start", comp, np.int64)
    copy_array(src_hap, dst_hap, "actual_length", "actual_length", comp, np.int64)

    seg_parent = readset_root.require_group("segments")
    seg = seg_parent.create_group(hap)
    copy_array(src_hap, seg, "read_index", "read_id", comp, np.int32)
    copy_array(src_hap, seg, "read_start", "hap_start", comp, np.int64)
    copy_array(src_hap, seg, "read_end", "hap_end", comp, np.int64)
    copy_array(src_hap, seg, "query_start", "query_start", comp, np.int32)
    copy_array(src_hap, seg, "query_end", "query_end", comp, np.int32)
    copy_array(src_hap, seg, "read_strand", "strand", comp, np.int8)
    copy_array(src_hap, seg, "window_read_offsets", "window_read_offsets", comp, np.int64)
    return {
        "hap": hap,
        "window_count": int(src_hap["sequence"].shape[0]),
        "total_bases": int(dst_hap["lengths"][:].sum()),
        "invalid_base_policy": invalid_policy,
        "invalid_base_count": int(seq_summary["invalid_base_count"]),
        "invalid_base_symbols": seq_summary["invalid_base_symbols"],
        "raw_start": src_hap["raw_start"][:].astype(np.int64).tolist(),
        "actual_length": src_hap["actual_length"][:].astype(np.int64).tolist(),
        "sequence_offsets": dst_hap["offsets"][:].astype(np.int64).tolist(),
        **copied_attrs,
    }


def selected_hap_attrs(src_hap: zarr.Group) -> dict[str, Any]:
    """Copy coordinate/window provenance from raw shared-HAP into compact output."""
    out: dict[str, Any] = {}
    for key in HAP_ATTR_COPY_KEYS:
        if key in src_hap.attrs:
            out[key] = src_hap.attrs[key]
    return out


def convert_reads_group(src_root: zarr.Group, src_reads: zarr.Group, readset_root: zarr.Group, dst_reads: zarr.Group, comp: Blosc, chunk_records: int) -> dict[str, Any]:
    dst_reads.attrs.update({
        "encoding": "4bit_acgtnp",
        "symbols": format_symbols(READ_ALLOWED),
        "kind": "read_sequence",
        "invalid_base_policy": "fail",
    })
    seq_summary = pack_object_sequence_array(src_reads["read_sequence"], dst_reads, READ_ALLOWED, comp, chunk_records, "fail")
    copy_array(src_reads, dst_reads, "read_haplotype", "hp", comp, np.uint8)
    copy_array(src_reads, dst_reads, "read_phase_set", "ps", comp, np.int64)
    has_cigar = "cigar" in src_root and "ops_packed" in src_root["cigar"]
    if "read_reference_start" in src_reads:
        copy_array(src_reads, dst_reads, "read_reference_start", "ref_start", comp, np.int64)
    if "read_reference_end" in src_reads:
        copy_array(src_reads, dst_reads, "read_reference_end", "ref_end", comp, np.int64)
    if "read_segment_type" in src_reads:
        copy_array(src_reads, dst_reads, "read_segment_type", "segment_type", comp, np.int8)
    if "read_cigar_offset" in src_reads:
        copy_array(src_reads, dst_reads, "read_cigar_offset", "cigar_offset", comp, np.int64)
    if "read_cigar_len" in src_reads:
        copy_array(src_reads, dst_reads, "read_cigar_len", "cigar_len", comp, np.int32)
    if has_cigar:
        cigar_group = readset_root.create_group("cigar")
        copy_array(src_root["cigar"], cigar_group, "ops_packed", "ops_packed", comp, np.uint32)
    return {
        "read_count": int(dst_reads["lengths"].shape[0]),
        "total_bases": int(dst_reads["lengths"][:].sum()),
        "packed_bytes": int(dst_reads["bases_packed"].shape[0]),
        "invalid_base_count": int(seq_summary["invalid_base_count"]),
        "invalid_base_symbols": seq_summary["invalid_base_symbols"],
        "has_cigar": bool(has_cigar),
        "cigar_op_count": int(src_root["cigar"]["ops_packed"].shape[0]) if has_cigar else 0,
    }


def pack_object_sequence_array(src_array: zarr.Array, dst_group: zarr.Group, allowed: frozenset[str], comp: Blosc, chunk_records: int, invalid_policy: str) -> dict[str, Any]:
    if invalid_policy not in {"fail", "to_n"}:
        raise ValueError(f"unknown invalid base policy: {invalid_policy}")
    bases = dst_group.empty("bases_packed", shape=(0,), dtype=np.uint8, chunks=(1 << 20,), compressor=comp)
    offsets = dst_group.empty("offsets", shape=(1,), dtype=np.int64, chunks=(8192,), compressor=comp)
    lengths = dst_group.empty("lengths", shape=(0,), dtype=np.int32, chunks=(8192,), compressor=comp)
    offsets[0] = 0

    base_offset = 0
    record_count = 0
    total_bases = 0
    invalid_counts: dict[str, int] = {}
    for lo in range(0, src_array.shape[0], chunk_records):
        hi = min(src_array.shape[0], lo + chunk_records)
        seqs = [str(value).upper() for value in src_array[lo:hi]]
        packed, local_offsets, local_lengths, local_total, local_invalid = pack_4bit_acgtnp_chunk(seqs, allowed, lo, invalid_policy)
        append_1d(bases, packed)
        append_1d(lengths, local_lengths)
        append_1d(offsets, local_offsets[1:] + base_offset)
        base_offset += int(packed.shape[0])
        record_count += len(seqs)
        total_bases += local_total
        add_counts(invalid_counts, local_invalid)
    dst_group.attrs.update({
        "record_count": record_count,
        "total_bases": total_bases,
        "packed_bytes": base_offset,
        "bases_per_byte": 2,
        "invalid_base_policy": invalid_policy,
        "invalid_base_count": int(sum(invalid_counts.values())),
        "invalid_base_symbols": dict(sorted(invalid_counts.items())),
    })
    return {
        "record_count": record_count,
        "total_bases": total_bases,
        "packed_bytes": base_offset,
        "invalid_base_count": int(sum(invalid_counts.values())),
        "invalid_base_symbols": dict(sorted(invalid_counts.items())),
    }


def pack_4bit_acgtnp_chunk(sequences, allowed, global_start_index, invalid_policy):
    from .compact.packing import pack_sequences
    return pack_sequences(sequences, allowed, global_start_index, invalid_policy,
                          ASCII_TO_CODE, format_symbols(allowed), _pack_4bit_acgtnp_chunk_legacy)


def _pack_4bit_acgtnp_chunk_legacy(sequences: list[str], allowed: frozenset[str], global_start_index: int, invalid_policy: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, dict[str, int]]:
    chunks: list[np.ndarray] = []
    offsets = np.zeros(len(sequences) + 1, dtype=np.int64)
    lengths = np.zeros(len(sequences), dtype=np.int32)
    byte_offset = 0
    total_bases = 0
    invalid_counts: dict[str, int] = {}
    allowed_codes = {int(ASCII_TO_CODE[ord(base)]) for base in allowed}
    n_code = int(ASCII_TO_CODE[ord("N")])
    for i, sequence in enumerate(sequences):
        raw = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
        codes = ASCII_TO_CODE[raw].copy()
        invalid = (codes == 255) | ~np.isin(codes, list(allowed_codes))
        if bool(invalid.any()):
            if invalid_policy == "fail":
                pos = int(np.flatnonzero(invalid)[0])
                base = chr(int(raw[pos]))
                raise ValueError(
                    f"invalid base for 4bit_acgtnp: base={base!r} record_index={global_start_index + i} "
                    f"position={pos} allowed={format_symbols(allowed)}"
                )
            values, counts = np.unique(raw[invalid], return_counts=True)
            for value, count in zip(values, counts, strict=False):
                key = chr(int(value)).upper()
                invalid_counts[key] = invalid_counts.get(key, 0) + int(count)
            codes[invalid] = n_code
        lengths[i] = raw.shape[0]
        offsets[i] = byte_offset
        if codes.shape[0] % 2:
            codes = np.concatenate([codes, np.zeros(1, dtype=np.uint8)])
        packed = ((codes[0::2] << 4) | codes[1::2]).astype(np.uint8, copy=False)
        chunks.append(packed)
        byte_offset += int(packed.shape[0])
        total_bases += int(lengths[i])
    offsets[len(sequences)] = byte_offset
    if chunks:
        return np.concatenate(chunks), offsets, lengths, total_bases, invalid_counts
    return np.asarray([], dtype=np.uint8), offsets, lengths, total_bases, invalid_counts


def add_counts(target: dict[str, int], source: dict[str, int]) -> None:
    for key, value in source.items():
        target[key] = target.get(key, 0) + int(value)


def merge_count_dicts(items: Iterable[dict[str, int]]) -> dict[str, int]:
    merged: dict[str, int] = {}
    for item in items:
        add_counts(merged, item)
    return dict(sorted(merged.items()))

def append_1d(array: zarr.Array, values: np.ndarray) -> None:
    if values.size == 0:
        return
    old = array.shape[0]
    array.resize(old + values.shape[0])
    array[old:old + values.shape[0]] = values


def copy_array(src_group: zarr.Group, dst_group: zarr.Group, src_name: str, dst_name: str, comp: Blosc, dtype: np.dtype) -> None:
    data = src_group[src_name][:].astype(dtype, copy=False)
    chunks = src_group[src_name].chunks
    if not chunks:
        chunks = (min(max(1, data.shape[0]), 8192),)
    dst_group.create_dataset(dst_name, data=data, chunks=chunks, compressor=comp, overwrite=True)


def make_window_rows(contig: str, hap: str, summary: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for i in range(summary["window_count"]):
        start = int(summary["sequence_offsets"][i])
        end = int(summary["sequence_offsets"][i + 1])
        rows.append({
            "contig": contig,
            "haplotype": hap,
            "window_index": i,
            "raw_start": int(summary["raw_start"][i]),
            "actual_length": int(summary["actual_length"][i]),
            "sequence_offset": start,
            "sequence_packed_bytes": end - start,
        })
    return rows


def write_windows_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    cols = ["contig", "haplotype", "window_index", "raw_start", "actual_length", "sequence_offset", "sequence_packed_bytes"]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(cols) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[col]) for col in cols) + "\n")


def write_contig_metadata(target_root: Path, config: LegacyConvertConfig, contig: str, source: Path, windows: list[dict[str, Any]], hap_summaries: dict[str, Any], read_summary: dict[str, Any]) -> None:
    payload = {
        "schema": "mfm_integrated_legacy_shared_hap_4bit_acgtnp_contig_v0",
        "sample_id": config.sample_id,
        "source_sample_id": config.source_sample_id,
        "dataset_role": config.dataset_role,
        "read_type": config.read_type,
        "project_id": config.project_id,
        "layout_name": config.layout_name,
        "contig": contig,
        "source_legacy_zarr": str(source),
        "window_rows": len(windows),
        "hap_summaries": hap_summaries,
        "read_summary": read_summary,
        "encoding": "4bit_acgtnp",
        "hap_invalid_policy": config.hap_invalid_policy,
    }
    write_json(target_root / "contig_metadata.json", payload)


def write_metadata(config: LegacyConvertConfig, summary: dict[str, Any]) -> None:
    metadata_dir = config.output_sample_root / "metadata"
    logs_dir = config.output_sample_root / "logs"
    tracks_dir = config.output_sample_root / "tracks" / "rna_bigwig"
    write_json(metadata_dir / "sample.json", {
        "sample_id": config.sample_id,
        "source_sample_id": config.source_sample_id,
        "dataset_role": config.dataset_role,
        "read_type": config.read_type,
        "project_id": config.project_id,
        "layout_name": config.layout_name,
        "schema": "mfm_integrated_sample_v0",
        "hap_invalid_policy": config.hap_invalid_policy,
    })
    write_json(metadata_dir / "preprocess.json", summary)
    write_tsv(metadata_dir / "input_manifest.tsv", [
        *[{"kind": "legacy_zarr", "path": str(config.legacy_root)}],
        *[{"kind": "source_wgs", "path": str(path)} for path in config.source_wgs],
        *[{"kind": "source_vcf", "path": str(path)} for path in config.source_vcfs],
        *[{"kind": "rna_bigwig", "path": str(path)} for path in config.rna_bigwigs],
    ], ["kind", "path"])
    outputs = []
    for contig_summary in summary["contigs"]:
        outputs.append({"kind": "compact_contig", "path": contig_summary.get("target", ""), "status": contig_summary.get("status", "")})
    write_tsv(metadata_dir / "output_manifest.tsv", outputs, ["kind", "path", "status"])
    write_json(logs_dir / "timing.json", {"timing": summary["timing"]})
    write_commands(metadata_dir / "commands.sh", config)
    if config.rna_bigwigs:
        write_rna_track_links(tracks_dir, config.rna_bigwigs)


def write_rna_track_links(tracks_dir: Path, bigwigs: list[Path]) -> None:
    tracks_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for index, source in enumerate(bigwigs):
        assay_id = source.stem.replace(".", "_")
        target_dir = tracks_dir / assay_id
        target_dir.mkdir(parents=True, exist_ok=True)
        target = target_dir / source.name
        if target.exists() or target.is_symlink():
            target.unlink()
        os.symlink(source, target)
        rows.append({"track_id": assay_id, "source_path": str(source), "linked_path": str(target)})
    write_tsv(tracks_dir / "track_manifest.tsv", rows, ["track_id", "source_path", "linked_path"])


def write_commands(path: Path, config: LegacyConvertConfig) -> None:
    contigs = "" if config.contigs is None else " --contigs " + ",".join(config.contigs)
    path.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "MFM-DNARNAFM/pre_training/preprocess_integrate_WGS_RNA/scripts/run_integrated_preprocess.sh "
        "convert-legacy-shared-hap "
        f"--sample-id {config.sample_id} "
        f"--source-sample-id {config.source_sample_id} "
        f"--legacy-root {config.legacy_root} "
        f"--output-sample-root {config.output_sample_root} "
        f"--dataset-role {config.dataset_role} "
        f"--read-type {config.read_type} "
        f"--project-id {config.project_id} "
        f"--hap-invalid-policy {config.hap_invalid_policy} "
        f"--layout-name {config.layout_name}"
        f"{contigs}\n",
        encoding="utf-8",
    )


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, Any]], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(cols) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(col, "")) for col in cols) + "\n")


def format_symbols(symbols: frozenset[str]) -> str:
    ordered = [base for base in SYMBOL_ORDER if base in symbols]
    extras = sorted(base for base in symbols if base not in SYMBOL_ORDER)
    return ",".join(ordered + extras)


def directory_size(path: Path) -> int:
    total = 0
    if path.is_file():
        return path.stat().st_size
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += (Path(root) / name).stat().st_size
            except FileNotFoundError:
                pass
    return total


def human_bytes(size: int) -> str:
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.2f} {unit}"
        value /= 1024
    return f"{value:.2f} TiB"
