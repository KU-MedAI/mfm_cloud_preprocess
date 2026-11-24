from __future__ import annotations

from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc, VLenUTF8


def write_store(path: Path, windows, reads: list[str], alignments: list[dict], associations: list[list[int]], attrs: dict) -> None:
    root = zarr.open_group(str(path), mode="w")
    root.attrs.update(attrs)
    text_codec = VLenUTF8()
    comp = Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE)
    wg = root.create_group("windows")
    wg.array("sequence", np.asarray([w.sequence for w in windows], dtype=object), object_codec=text_codec)
    for name, values in {
        "raw_start": [w.raw_start for w in windows], "actual_length": [w.actual_length for w in windows]
    }.items():
        wg.array(name, np.asarray(values, dtype=np.int64), compressor=comp)
    rg = root.create_group("reads")
    rg.array("sequence", np.asarray(reads, dtype=object), object_codec=text_codec)
    ag = root.create_group("alignments")
    fields = {
        "read_index": np.int64, "raw_start": np.int64, "raw_end": np.int64,
        "reference_start": np.int64, "reference_end": np.int64,
        "query_start": np.int32, "query_end": np.int32, "strand": np.int8,
        "segment_type": np.int8,
    }
    for name, dtype in fields.items():
        ag.array(name, np.asarray([a[name] for a in alignments], dtype=dtype), compressor=comp)
    ag.array("cigar", np.asarray([a["cigar"] for a in alignments], dtype=object), object_codec=text_codec)
    flat = [a for group in associations for a in group]
    offsets = np.cumsum([0] + [len(group) for group in associations], dtype=np.int64)
    xg = root.create_group("window_alignments")
    xg.array("alignment_indices", np.asarray(flat, dtype=np.int64), compressor=comp)
    xg.array("offsets", offsets, compressor=comp)


def overlapping_alignment_indices(root, window_index: int, crop_start: int, crop_size: int = 8192) -> np.ndarray:
    offsets = root["window_alignments/offsets"][:]
    idx = root["window_alignments/alignment_indices"][offsets[window_index]:offsets[window_index + 1]]
    window_start = int(root["windows/raw_start"][window_index])
    start, end = window_start + crop_start, window_start + crop_start + crop_size
    return idx[(root["alignments/raw_start"][idx] < end) & (root["alignments/raw_end"][idx] > start)]
