"""기존 read 필터를 공유하는 CIGAR 구간 집계와 별도 aligned-pairs 검증기임."""
from collections import OrderedDict
import os
import numpy as np
from pre_training.preprocess_raw_readset.filters import rejection_reason

BLOCK = 262144
TABLE = np.full(256, 255, dtype=np.uint8)
TABLE[ord('N')] = TABLE[ord('n')] = 4
for i, base in enumerate(b'ACGT'):
    TABLE[base] = TABLE[base + 32] = i


class Counts:
    def __init__(self, start, end):
        self.start, self.end = start, end
        n = end - start
        self.bases = np.zeros((5, n), dtype=np.uint32)
        self.mapq = np.zeros(n, dtype=np.uint64)
        self.baseq = np.zeros(n, dtype=np.uint64)
        self.baseq_n = np.zeros(n, dtype=np.uint32)
        self.forward = np.zeros(n, dtype=np.uint32)
        self.deletion = np.zeros(n, dtype=np.uint32)
        self.insertion = np.zeros(n, dtype=np.uint32)
        self.deletion_event = np.zeros(n, dtype=np.uint32)
        self.max_deletion = np.zeros(n, dtype=np.uint32)
        self.supplementary = np.zeros(n, dtype=np.uint32)

    def observations(self, read, refs, queries, codes, qualities):
        index = refs - self.start
        if read.is_supplementary:
            self.supplementary[index] += 1
            return
        values = codes[queries]
        if np.any(values > 4):
            raise ValueError('관측 read에 지원하지 않는 염기 문자가 있음')
        for base in range(5):
            self.bases[base, index] += values == base
        self.mapq[index] += int(read.mapping_quality)
        self.forward[index] += not read.is_reverse
        if qualities is not None:
            q = qualities[queries]
            valid = q != 255
            self.baseq[index] += np.where(valid, q, 0).astype(np.uint64)
            self.baseq_n[index] += valid

    def event(self, anchor, length, insertion):
        if not self.start <= anchor < self.end:
            return
        i = anchor - self.start
        if insertion:
            self.insertion[i] += 1
        else:
            self.deletion_event[i] += 1
            self.max_deletion[i] = max(int(self.max_deletion[i]), length)

    def finish(self):
        total = self.bases.sum(axis=0, dtype=np.uint64)
        if np.any(total > np.iinfo(np.uint32).max):
            raise OverflowError('uint32 depth 범위 초과')
        depth = total.astype(np.uint32)
        denominator = total + self.deletion.astype(np.uint64)
        features = np.zeros((len(depth), 12), dtype=np.float32)
        features[:, 0] = np.log1p(total.astype(np.float64))
        def ratio(values, divisor):
            return np.divide(values, divisor, out=np.zeros(len(depth), dtype=np.float64), where=divisor > 0)
        features[:, 1] = ratio(self.mapq, total)
        features[:, 2] = ratio(self.baseq, self.baseq_n)
        features[:, 3] = ratio(self.forward, total)
        features[:, 4:9] = np.stack([ratio(v, total) for v in self.bases], axis=1)
        features[:, 9] = ratio(self.insertion, denominator)
        features[:, 10] = ratio(self.deletion_event, denominator)
        features[:, 11] = self.max_deletion
        if not np.isfinite(features).all():
            raise ValueError('pileup feature의 NaN/inf')
        return dict(features=features, depth=depth, baseq_depth=self.baseq_n,
                    supplementary_depth=self.supplementary, max_deletion_length=self.max_deletion)


def calculate(bam, contig, start, end, config, *, oracle=False):
    counts = Counts(start, end)
    for read in bam.fetch(contig, start, end):
        if rejection_reason(read, config):
            continue
        sequence = np.frombuffer(read.query_sequence.encode('ascii'), dtype=np.uint8)
        codes = TABLE[sequence]
        qualities = np.asarray(read.query_qualities, dtype=np.uint8) if read.query_qualities is not None else None
        if oracle:
            # HTSlib의 독립 query↔reference 대응을 사용하여 production CIGAR walker와 대조함.
            pairs = read.get_aligned_pairs(with_cigar=True)
            if not pairs:
                continue
            q = np.fromiter((a if a is not None else -1 for a, _, _ in pairs), dtype=np.int64)
            r = np.fromiter((b if b is not None else -1 for _, b, _ in pairs), dtype=np.int64)
            ops = np.fromiter((int(c) for _, _, c in pairs), dtype=np.int8)
            # HTSlib aligned-pairs의 P query 증가를 SAM의 query 비소모 규칙으로 보정함.
            q = np.where(q >= 0, q - np.cumsum(ops == 6), -1)
            chosen = (r >= start) & (r < end) & (q >= 0) & np.isin(ops, (0, 7, 8))
            counts.observations(read, r[chosen], q[chosen], codes, qualities)
            if read.is_supplementary:
                continue
            deleted = (r >= start) & (r < end) & (ops == 2)
            counts.deletion[r[deleted] - start] += 1
            boundaries = np.r_[0, np.flatnonzero(ops[1:] != ops[:-1]) + 1, len(ops)]
            previous = np.r_[read.reference_start-1, np.maximum.accumulate(r)]
            starts, lengths = boundaries[:-1], np.diff(boundaries)
            for op in (1, 2):
                chosen = ops[starts] == op
                anchors = previous[starts[chosen]] if op == 1 else r[starts[chosen]]-1
                size = lengths[chosen]
                valid = (anchors >= max(start, read.reference_start)) & (anchors < end)
                index = anchors[valid]-start
                np.add.at(counts.insertion if op == 1 else counts.deletion_event, index, 1)
                if op == 2:
                    np.maximum.at(counts.max_deletion, index, size[valid].astype(np.uint32))
            continue
        ref, query = read.reference_start, 0
        matches = []
        for op, length in read.cigartuples or []:
            if op in (0, 7, 8):
                a, b = max(ref, start), min(ref + length, end)
                if a < b:
                    matches.append((a, b-a, query-ref))
                ref += length
                query += length
            elif op == 1:
                if not read.is_supplementary and ref > read.reference_start:
                    counts.event(ref-1, length, True)
                query += length
            elif op in (2, 3):
                if op == 2 and not read.is_supplementary:
                    a, b = max(ref, start), min(ref + length, end)
                    if a < b:
                        counts.deletion[a-start:b-start] += 1
                    if ref > read.reference_start:
                        counts.event(ref-1, length, False)
                ref += length
            elif op == 4:
                query += length
            elif op not in (5, 6):
                raise ValueError('지원하지 않는 CIGAR op: ' + str(op))
        if matches:
            # 짧은 CIGAR run마다 집계하지 않고, block 안의 한 read 관측을 한 번에 집계함.
            starts, lengths, shifts = np.asarray(matches, dtype=np.int64).T
            offsets = np.r_[0, np.cumsum(lengths)[:-1]]
            refs = np.arange(int(lengths.sum())) + np.repeat(starts-offsets, lengths)
            queries = refs + np.repeat(shifts, lengths)
            counts.observations(read, refs, queries, codes, qualities)
    return counts.finish()


class ReferenceCache:
    def __init__(self, bam, contig, config, *, oracle=False, block=BLOCK, cache_blocks=12, validation_backend=None):
        self.bam, self.contig, self.config = bam, contig, config
        self.oracle, self.block, self.limit = oracle, block, cache_blocks
        self.length = bam.get_reference_length(contig)
        self.cache = OrderedDict()
        self.native = None
        self.blocks_computed = self.cache_hits = 0
        self.backend = 'production_cigar'
        if oracle:
            requested = validation_backend or os.environ.get('MFM_PILEUP_VALIDATION_BACKEND','auto')
            if requested not in {'auto','native','aligned_pairs'}:
                raise ValueError('MFM_PILEUP_VALIDATION_BACKEND는 auto/native/aligned_pairs임')
            self.backend = 'aligned_pairs'
            if requested != 'aligned_pairs':
                from ..record_stats_parallel import NativeUnavailable
                from .native_oracle import NativeOracle
                try:
                    self.native = NativeOracle(bam,contig,config)
                    self.backend = 'native_htslib_pileup'
                except NativeUnavailable as error:
                    if requested == 'native':
                        raise
                    print('[PILEUP-VALIDATION] aligned_pairs fallback: '+str(error),flush=True)
        else:
            requested = os.environ.get('MFM_PILEUP_PRODUCTION_BACKEND','auto')
            if requested not in {'auto','native','python'}:
                raise ValueError('MFM_PILEUP_PRODUCTION_BACKEND는 auto/native/python임')
            if requested != 'python':
                from ..record_stats_parallel import NativeUnavailable
                from .native_oracle import NativeProducer
                try:
                    self.native = NativeProducer(bam,contig,config)
                    self.backend = 'native_cigar'
                except NativeUnavailable as error:
                    if requested == 'native':raise
                    print('[PILEUP-PRODUCTION] Python CIGAR fallback: '+str(error),flush=True)

    def close(self):
        if self.native is not None:
            self.native.close()
        self.cache.clear()

    def gather(self, positions, names=None):
        names = names or ('features', 'depth', 'baseq_depth', 'supplementary_depth', 'max_deletion_length')
        shape = len(positions)
        out = {n: np.zeros((shape, 12) if n == 'features' else shape, dtype=np.float32 if n == 'features' else np.uint32) for n in names}
        valid = (positions >= 0) & (positions < self.length)
        blocks = np.unique(positions[valid] // self.block)
        for block in blocks:
            key = int(block)
            value = self.cache.pop(key, None)
            if value is None:
                start = key * self.block
                end = min(start+self.block,self.length)
                value = self.native.calculate(start,end) if self.native is not None else calculate(self.bam,self.contig,start,end,self.config,oracle=self.oracle)
                self.blocks_computed += 1
            else:
                self.cache_hits += 1
            self.cache[key] = value
            while len(self.cache) > self.limit:
                self.cache.popitem(last=False)
            chosen = valid & (positions // self.block == key)
            index = positions[chosen] - key*self.block
            for name in names:
                out[name][chosen] = value[name][index]
        return out


def projected(cache, projection, start, actual, padded):
    positions, anchors, kind = projection.window(start, actual, padded)
    result = cache.gather(positions)
    inserted = kind == 2
    if inserted.any():
        # 삽입 염기의 관측은 만들지 않고 anchor의 generic 삽입 사건 비율만 전달함.
        anchor = cache.gather(anchors[inserted], ('features',))
        result['features'][inserted, 9] = anchor['features'][:, 9]
    result['mapping_kind'] = kind
    return result
