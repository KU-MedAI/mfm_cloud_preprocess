from __future__ import annotations
from pathlib import Path
import numpy as np
import zarr
from numcodecs import Blosc, VLenUTF8

def _records_from_pairs(windows, reads, alignments, associations):
    """Build compact window records for compatibility with legacy callers."""
    records = []
    for wi, indices in enumerate(associations):
        window_records = []
        for ai in indices:
            alignment = alignments[ai]
            full = reads[alignment["read_index"]]
            pairs = alignment.get("projected_pairs", [])
            w0 = windows[wi].raw_start
            w1 = w0 + windows[wi].actual_length
            selected = [(q, raw) for q, raw in pairs if w0 <= raw < w1]
            if not selected:
                continue
            sequence = "".join(full[q] for q, _ in selected)
            window_records.append({
                "sequence": sequence,
                "start": min(raw for _, raw in selected) - w0,
                "end": max(raw for _, raw in selected) + 1 - w0,
                "strand": int(alignment["strand"]),
                "alignment_index": int(ai),
            })
        records.append(window_records)
    return records


def write_model_store(
    path: Path,
    windows,
    reads: list[str],
    alignments: list[dict],
    associations: list[list[int]],
    attrs: dict,
    *,
    model_records: list[list[dict]] | None = None,
) -> None:
    root = zarr.open_group(str(path), mode="w")
    root.attrs.update({"schema": "model_input_v3", "base_encoding": "literal_string", **attrs})
    comp = Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE)
    root.array("sequence", np.asarray([w.sequence for w in windows], dtype=object), object_codec=VLenUTF8())
    root.array("raw_start", np.asarray([w.raw_start for w in windows], dtype=np.int64), compressor=comp)
    root.array("actual_length", np.asarray([w.actual_length for w in windows], dtype=np.int64), compressor=comp)
    if model_records is None:
        model_records = _records_from_pairs(windows, reads, alignments, associations)
    if len(model_records) != len(windows):
        raise ValueError("model_records must contain one list per window")
    packed=[]; seq_offsets=[0]; starts=[]; ends=[]; lengths=[]; strands=[]; alignment_indices=[]; window_offsets=[0]
    for records in model_records:
        for record in records:
            seq = str(record["sequence"])
            if not seq:
                continue
            packed.append(seq)
            seq_offsets.append(seq_offsets[-1] + 1)
            starts.append(int(record["start"]))
            ends.append(int(record["end"]))
            lengths.append(len(seq))
            strands.append(int(record["strand"]))
            alignment_indices.append(int(record["alignment_index"]))
        window_offsets.append(len(starts))
    flat=np.asarray(packed, dtype=object) if packed else np.empty(0,dtype=object)
    root.array("read_sequence",flat,object_codec=VLenUTF8()); root.array("read_sequence_offsets",np.asarray(seq_offsets,dtype=np.int64),compressor=comp)
    root.array("read_start",np.asarray(starts,dtype=np.int64),compressor=comp); root.array("read_end",np.asarray(ends,dtype=np.int64),compressor=comp)
    root.array("read_length",np.asarray(lengths,dtype=np.int32),compressor=comp); root.array("read_strand",np.asarray(strands,dtype=np.int8),compressor=comp)
    root.array("alignment_index",np.asarray(alignment_indices,dtype=np.int64),compressor=comp)
    root.array("window_read_offsets",np.asarray(window_offsets,dtype=np.int64),compressor=comp)

def load_window(root, window_index: int):
    offsets=root["window_read_offsets"][:]; lo,hi=int(offsets[window_index]),int(offsets[window_index+1]); so=root["read_sequence_offsets"][lo:hi+1]; raw=root["read_sequence"]
    seqs=[]
    for i in range(hi-lo):
        value=raw[int(so[i])]
        seqs.append(value.item() if hasattr(value, "item") else str(value))
    return root["sequence"][window_index], {"sequence":seqs,"start":root["read_start"][lo:hi],"end":root["read_end"][lo:hi],"length":root["read_length"][lo:hi],"strand":root["read_strand"][lo:hi]}
