from __future__ import annotations

import numpy as np


def pack_cigar(cigartuples: list[tuple[int, int]] | None) -> list[int]:
    if not cigartuples:
        return []
    packed: list[int] = []
    for op_code, length in cigartuples:
        if length < 0:
            raise ValueError(f"negative CIGAR length: {length}")
        if length >= (1 << 28):
            raise ValueError(f"CIGAR length too large for uint32 packing: {length}")
        packed.append((int(length) << 4) | int(op_code))
    return packed


def to_uint32(values: list[int]) -> np.ndarray:
    return np.asarray(values, dtype=np.uint32)
