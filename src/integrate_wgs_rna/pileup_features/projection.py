"""검증된 consensus chain을 bounded HAP 구간에만 펼침. 염색체 전체 배열은 만들지 않음."""
from pathlib import Path
import numpy as np
from pre_training.preprocess_raw_readset.coordinates import load_chain_alignment


class Projection:
    def __init__(self, path, contig, reference_length, hap_length):
        self.chain = load_chain_alignment(Path(path), query_name=contig, target_name=contig)
        h = self.chain.header
        if (h['target_size'], h['query_size']) != (reference_length, hap_length):
            raise ValueError('chain과 reference/HAP 길이가 다름')
        if (h['target_start'], h['query_start'], h['target_end'], h['query_end']) != (0, 0, reference_length, hap_length):
            raise ValueError('완전한 plus/plus consensus chain이 필요함')
        self.blocks = self.chain.blocks
        self.ends = np.array([q + n for q, _, n in self.blocks], dtype=np.int64)
        self.insertions = []
        q0 = t0 = 0
        for q, t, n in self.blocks:
            if n <= 0 or q < q0 or t < t0:
                raise ValueError('chain block 순서/길이 오류')
            if q > q0:
                self.insertions.append((q0, q, t0 - 1))
            q0, t0 = q + n, t + n
        if q0 < hap_length:
            self.insertions.append((q0, hap_length, t0 - 1))
        self.insertion_ends = np.array([b for _, b, _ in self.insertions], dtype=np.int64)

    def window(self, start, actual, padded):
        if start < 0 or not 0 < actual <= padded or start + actual > self.chain.header['query_size']:
            raise ValueError('HAP window 범위 오류')
        positions = np.full(padded, -1, dtype=np.int64)
        anchors = np.full(padded, -1, dtype=np.int64)
        # 0=padding, 1=reference 대응, 2=삽입. 깊이0과 비대응을 구분함.
        kind = np.zeros(padded, dtype=np.uint8)
        stop = start + actual
        first = int(np.searchsorted(self.ends, start, side='right'))
        for index in range(first, len(self.blocks)):
            q, t, n = self.blocks[index]
            if q >= stop:
                break
            a, b = max(start, q), min(stop, q + n)
            if a < b:
                positions[a-start:b-start] = np.arange(t + a-q, t + b-q, dtype=np.int64)
                kind[a-start:b-start] = 1
        first = int(np.searchsorted(self.insertion_ends, start, side='right'))
        for index in range(first, len(self.insertions)):
            a, b, anchor = self.insertions[index]
            if a >= stop:
                break
            a, b = max(a, start), min(b, stop)
            if a < b:
                anchors[a-start:b-start] = anchor
                kind[a-start:b-start] = 2
        if np.any(kind[:actual] == 0):
            raise ValueError('chain이 실제 HAP 구간을 설명하지 못함')
        return positions, anchors, kind
