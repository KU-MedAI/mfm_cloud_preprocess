from __future__ import annotations

from dataclasses import dataclass

import numpy as np


ACGT_TO_2BIT = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
}
TWO_BIT_TO_BASE = "ACGT"

IUPAC_TO_4BIT = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "N": 4,
    "R": 5,
    "Y": 6,
    "S": 7,
    "W": 8,
    "K": 9,
    "M": 10,
    "B": 11,
    "D": 12,
    "H": 13,
    "V": 14,
    "P": 15,
}
FOUR_BIT_TO_BASE = {
    code: base for base, code in IUPAC_TO_4BIT.items()
}

ACGTNP_TO_4BIT = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "N": 4,
    "P": 5,
}
FOUR_BIT_ACGTNP_TO_BASE = {
    code: base for base, code in ACGTNP_TO_4BIT.items()
}
ACGTNP_SYMBOLS = frozenset(ACGTNP_TO_4BIT)
ACGTNP_SYMBOL_ORDER = ("A", "C", "G", "T", "N", "P")


@dataclass(frozen=True)
class PackedReads:
    bases: np.ndarray
    offsets: np.ndarray
    lengths: np.ndarray
    sparse_offsets: np.ndarray | None
    sparse_positions: np.ndarray | None
    sparse_codes: np.ndarray | None
    attrs: dict


def count_non_acgt(sequence: str) -> tuple[int, int]:
    non_acgt = 0
    n_count = 0
    for base in sequence:
        if base not in ACGT_TO_2BIT:
            non_acgt += 1
            if base == "N":
                n_count += 1
    return non_acgt, n_count


def pack_reads(
    sequences: list[str],
    encoding: str,
    *,
    allowed_symbols: frozenset[str] | set[str] | None = None,
) -> PackedReads:
    if encoding == "2bit_acgt":
        return _pack_2bit_acgt(sequences)
    if encoding == "4bit_acgtnp":
        return _pack_4bit_acgtnp(sequences, allowed_symbols=allowed_symbols)
    if encoding == "4bit_iupac":
        raise ValueError("4bit_iupac is legacy decode-only; use strict 4bit_acgtnp for new data")
    if encoding == "2bit_acgt_n_sparse":
        return _pack_2bit_acgt_sparse(sequences, allowed_symbols=allowed_symbols)
    raise ValueError(f"unknown read encoding: {encoding}")


def decode_read(packed: PackedReads, read_index: int, encoding: str) -> str:
    if encoding == "2bit_acgt":
        return _decode_2bit(packed, read_index)
    if encoding == "4bit_acgtnp":
        return _decode_4bit_acgtnp(packed, read_index)
    if encoding == "4bit_iupac":
        return _decode_4bit_iupac(packed, read_index)
    if encoding == "2bit_acgt_n_sparse":
        return _decode_2bit_sparse(packed, read_index)
    raise ValueError(f"unknown read encoding: {encoding}")


def _pack_2bit_acgt(sequences: list[str]) -> PackedReads:
    offsets = np.zeros(len(sequences) + 1, dtype=np.int64)
    lengths = np.zeros(len(sequences), dtype=np.int32)
    chunks: list[bytes] = []
    non_acgt_bases = 0
    n_bases = 0
    non_acgt_reads = 0
    byte_offset = 0
    for i, sequence in enumerate(sequences):
        sequence = sequence.upper()
        lengths[i] = len(sequence)
        offsets[i] = byte_offset
        local_non_acgt, local_n = count_non_acgt(sequence)
        non_acgt_bases += local_non_acgt
        n_bases += local_n
        if local_non_acgt:
            non_acgt_reads += 1
        packed = bytearray((len(sequence) + 3) // 4)
        for pos, base in enumerate(sequence):
            code = ACGT_TO_2BIT.get(base, 0)
            packed[pos // 4] |= code << ((3 - (pos % 4)) * 2)
        chunks.append(bytes(packed))
        byte_offset += len(packed)
    offsets[len(sequences)] = byte_offset
    return PackedReads(
        bases=np.frombuffer(b"".join(chunks), dtype=np.uint8).copy(),
        offsets=offsets,
        lengths=lengths,
        sparse_offsets=None,
        sparse_positions=None,
        sparse_codes=None,
        attrs={
            "encoding": "2bit_acgt",
            "lossless": non_acgt_bases == 0,
            "non_acgt_bases_replaced_as_A": int(non_acgt_bases),
            "n_bases_replaced_as_A": int(n_bases),
            "non_acgt_reads": int(non_acgt_reads),
            "symbols": "A,C,G,T",
        },
    )


def _pack_2bit_acgt_sparse(
    sequences: list[str],
    *,
    allowed_symbols: frozenset[str] | set[str] | None = None,
) -> PackedReads:
    allowed = ACGTNP_SYMBOLS if allowed_symbols is None else frozenset(allowed_symbols)
    offsets = np.zeros(len(sequences) + 1, dtype=np.int64)
    lengths = np.zeros(len(sequences), dtype=np.int32)
    sparse_offsets = np.zeros(len(sequences) + 1, dtype=np.int64)
    sparse_positions: list[int] = []
    sparse_codes: list[int] = []
    chunks: list[bytes] = []
    byte_offset = 0
    non_acgt_bases = 0
    n_bases = 0
    for i, sequence in enumerate(sequences):
        sequence = sequence.upper()
        lengths[i] = len(sequence)
        offsets[i] = byte_offset
        sparse_offsets[i] = len(sparse_positions)
        packed = bytearray((len(sequence) + 3) // 4)
        for pos, base in enumerate(sequence):
            validate_symbol(base, allowed, "2bit_acgt_n_sparse", i, pos)
            code = ACGT_TO_2BIT.get(base)
            if code is None:
                code = 0
                sparse_positions.append(pos)
                sparse_codes.append(IUPAC_TO_4BIT.get(base, IUPAC_TO_4BIT["N"]))
                non_acgt_bases += 1
                if base == "N":
                    n_bases += 1
            packed[pos // 4] |= code << ((3 - (pos % 4)) * 2)
        chunks.append(bytes(packed))
        byte_offset += len(packed)
    offsets[len(sequences)] = byte_offset
    sparse_offsets[len(sequences)] = len(sparse_positions)
    return PackedReads(
        bases=np.frombuffer(b"".join(chunks), dtype=np.uint8).copy(),
        offsets=offsets,
        lengths=lengths,
        sparse_offsets=sparse_offsets,
        sparse_positions=np.asarray(sparse_positions, dtype=np.int32),
        sparse_codes=np.asarray(sparse_codes, dtype=np.uint8),
        attrs={
            "encoding": "2bit_acgt_n_sparse",
            "lossless": True,
            "non_acgt_bases_sparse": int(non_acgt_bases),
            "n_bases_sparse": int(n_bases),
            "symbols": format_symbols(allowed),
        },
    )


def _pack_4bit_acgtnp(
    sequences: list[str],
    *,
    allowed_symbols: frozenset[str] | set[str] | None = None,
) -> PackedReads:
    allowed = ACGTNP_SYMBOLS if allowed_symbols is None else frozenset(allowed_symbols)
    offsets = np.zeros(len(sequences) + 1, dtype=np.int64)
    lengths = np.zeros(len(sequences), dtype=np.int32)
    chunks: list[bytes] = []
    byte_offset = 0
    for i, sequence in enumerate(sequences):
        sequence = sequence.upper()
        lengths[i] = len(sequence)
        offsets[i] = byte_offset
        packed = bytearray((len(sequence) + 1) // 2)
        for pos, base in enumerate(sequence):
            validate_symbol(base, allowed, "4bit_acgtnp", i, pos)
            code = ACGTNP_TO_4BIT[base]
            if pos % 2 == 0:
                packed[pos // 2] |= code << 4
            else:
                packed[pos // 2] |= code
        chunks.append(bytes(packed))
        byte_offset += len(packed)
    offsets[len(sequences)] = byte_offset
    return PackedReads(
        bases=np.frombuffer(b"".join(chunks), dtype=np.uint8).copy(),
        offsets=offsets,
        lengths=lengths,
        sparse_offsets=None,
        sparse_positions=None,
        sparse_codes=None,
        attrs={
            "encoding": "4bit_acgtnp",
            "lossless": True,
            "reserved_codes": "6,7,8,9,10,11,12,13,14,15",
            "symbols": format_symbols(allowed),
        },
    )


def _pack_4bit_iupac(sequences: list[str]) -> PackedReads:
    # Legacy helper kept only for historical benchmark/reference code. New data must use 4bit_acgtnp.
    offsets = np.zeros(len(sequences) + 1, dtype=np.int64)
    lengths = np.zeros(len(sequences), dtype=np.int32)
    chunks: list[bytes] = []
    byte_offset = 0
    unknown_bases = 0
    for i, sequence in enumerate(sequences):
        sequence = sequence.upper()
        lengths[i] = len(sequence)
        offsets[i] = byte_offset
        packed = bytearray((len(sequence) + 1) // 2)
        for pos, base in enumerate(sequence):
            code = IUPAC_TO_4BIT.get(base)
            if code is None:
                code = IUPAC_TO_4BIT["N"]
                unknown_bases += 1
            if pos % 2 == 0:
                packed[pos // 2] |= code << 4
            else:
                packed[pos // 2] |= code
        chunks.append(bytes(packed))
        byte_offset += len(packed)
    offsets[len(sequences)] = byte_offset
    return PackedReads(
        bases=np.frombuffer(b"".join(chunks), dtype=np.uint8).copy(),
        offsets=offsets,
        lengths=lengths,
        sparse_offsets=None,
        sparse_positions=None,
        sparse_codes=None,
        attrs={
            "encoding": "4bit_iupac",
            "lossless": unknown_bases == 0,
            "unknown_bases_encoded_as_N": int(unknown_bases),
            "symbols": "A,C,G,T,N,R,Y,S,W,K,M,B,D,H,V,P",
        },
    )


def format_symbols(symbols: frozenset[str]) -> str:
    ordered = [base for base in ACGTNP_SYMBOL_ORDER if base in symbols]
    extras = sorted(base for base in symbols if base not in ACGTNP_SYMBOL_ORDER)
    return ",".join(ordered + extras)


def validate_symbol(base: str, allowed: frozenset[str], encoding: str, read_index: int, pos: int) -> None:
    if base not in allowed:
        allowed_text = format_symbols(allowed)
        raise ValueError(
            f"invalid base for {encoding}: base={base!r} read_index={read_index} "
            f"position={pos} allowed={allowed_text}"
        )


def _decode_2bit(packed: PackedReads, read_index: int) -> str:
    start = int(packed.offsets[read_index])
    length = int(packed.lengths[read_index])
    data = packed.bases[start:int(packed.offsets[read_index + 1])]
    out = []
    for pos in range(length):
        byte = int(data[pos // 4])
        shift = (3 - (pos % 4)) * 2
        out.append(TWO_BIT_TO_BASE[(byte >> shift) & 0b11])
    return "".join(out)


def _decode_2bit_sparse(packed: PackedReads, read_index: int) -> str:
    out = list(_decode_2bit(packed, read_index))
    if packed.sparse_offsets is None or packed.sparse_positions is None or packed.sparse_codes is None:
        return "".join(out)
    start = int(packed.sparse_offsets[read_index])
    end = int(packed.sparse_offsets[read_index + 1])
    for pos, code in zip(packed.sparse_positions[start:end], packed.sparse_codes[start:end], strict=False):
        out[int(pos)] = FOUR_BIT_TO_BASE.get(int(code), "N")
    return "".join(out)


def _decode_4bit_acgtnp(packed: PackedReads, read_index: int) -> str:
    start = int(packed.offsets[read_index])
    length = int(packed.lengths[read_index])
    data = packed.bases[start:int(packed.offsets[read_index + 1])]
    out = []
    for pos in range(length):
        byte = int(data[pos // 2])
        code = (byte >> 4) if pos % 2 == 0 else (byte & 0x0F)
        if code not in FOUR_BIT_ACGTNP_TO_BASE:
            raise ValueError(f"invalid reserved 4bit_acgtnp code={code} read_index={read_index} position={pos}")
        out.append(FOUR_BIT_ACGTNP_TO_BASE[code])
    return "".join(out)


def _decode_4bit_iupac(packed: PackedReads, read_index: int) -> str:
    start = int(packed.offsets[read_index])
    length = int(packed.lengths[read_index])
    data = packed.bases[start:int(packed.offsets[read_index + 1])]
    out = []
    for pos in range(length):
        byte = int(data[pos // 2])
        code = (byte >> 4) if pos % 2 == 0 else (byte & 0x0F)
        out.append(FOUR_BIT_TO_BASE.get(code, "N"))
    return "".join(out)
