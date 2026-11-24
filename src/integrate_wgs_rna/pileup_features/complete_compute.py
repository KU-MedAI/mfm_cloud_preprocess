"""공통/HP 관측을 한 번 읽기로 계산하고 HAP 좌표에 대응함."""
from collections import OrderedDict
import os
import numpy as np
import pysam
from pre_training.preprocess_raw_readset.filters import rejection_reason
from ..record_stats_parallel import NativeUnavailable
from .schema import ARRAYS, CHANNELS, COUNT_FIELDS, INSERTION_CHANNELS
from .variants import VariantLookup

INDEX = {name: i for i, name in enumerate(COUNT_FIELDS)}
CHANNEL_INDEX = {name: i for i, name in enumerate(CHANNELS)}
TABLE = np.full(256, 255, dtype=np.uint8)
for i, base in enumerate(b'ACGTN'):
    TABLE[base] = TABLE[base+32] = i


def calculate_counts(bam, contig, start, end, config):
    """Compiler가 없을 때도 같은 필터와 단일 fetch를 사용하는 경로임."""
    counts = np.zeros((3, len(COUNT_FIELDS), end-start), dtype=np.uint32)
    sums = np.zeros((3, 2, end-start), dtype=np.uint64)
    for read in bam.fetch(contig, start, end):
        if rejection_reason(read, config):
            continue
        tag = read.get_tag('HP') if read.has_tag('HP') else 0
        hp = tag if isinstance(tag, int) and tag in (1, 2) else 0
        groups = (0, hp) if hp else (0,)
        seq = TABLE[np.frombuffer(read.query_sequence.encode('ascii'), dtype=np.uint8)]
        qual = np.asarray(read.query_qualities, dtype=np.uint8) if read.query_qualities is not None else None
        ref, query = read.reference_start, 0
        # 같은 연산이 연속된 CIGAR는 한 사건으로 정규화함.
        cigar = []
        for op, size in read.cigartuples or []:
            if cigar and cigar[-1][0] == op:
                cigar[-1] = (op, cigar[-1][1]+size)
            else:
                cigar.append((op, size))
        for op, size in cigar:
            if size <= 0:
                raise ValueError('0 길이 CIGAR')
            if op in (0, 7, 8):
                a, b = max(ref, start), min(ref+size, end)
                if a < b:
                    p = np.arange(a-start, b-start)
                    q = np.arange(query+a-ref, query+b-ref)
                    if np.any(q < 0) or np.any(q >= len(seq)):
                        raise ValueError('CIGAR query 길이 불일치')
                    if read.is_supplementary:
                        for group in groups:
                            counts[group, INDEX['supplementary'], p] += 1
                    else:
                        observed = seq[q]
                        if np.any(observed > 4):
                            raise ValueError('지원하지 않는 관측 read 염기')
                        for group in groups:
                            for base in range(5):
                                counts[group, base, p] += observed == base
                            counts[group, INDEX['forward'], p] += not read.is_reverse
                        sums[0, 0, p] += int(read.mapping_quality)
                        if qual is not None:
                            valid = qual[q] != 255
                            counts[0, INDEX['baseq_n'], p] += valid
                            sums[0, 1, p] += np.where(valid, qual[q], 0).astype(np.uint64)
                ref += size
                query += size
            elif op in (1, 2):
                anchor = ref-1
                if not read.is_supplementary:
                    for group in groups:
                        if ref > read.reference_start and start <= anchor < end:
                            key = 'insertion' if op == 1 else 'deletion_event'
                            counts[group, INDEX[key], anchor-start] += 1
                            if op == 2:
                                j = INDEX['max_deletion']
                                counts[group, j, anchor-start] = max(int(counts[group, j, anchor-start]), size)
                        if op == 2:
                            a, b = max(ref, start), min(ref+size, end)
                            if a < b:
                                counts[group, INDEX['deletion'], a-start:b-start] += 1
                if op == 1:
                    query += size
                else:
                    ref += size
            elif op == 3:
                ref += size
            elif op == 4:
                query += size
            elif op not in (5, 6):
                raise ValueError('지원하지 않는 CIGAR')
        if query != len(seq):
            raise ValueError('CIGAR query 길이 불일치')
    return counts, sums


def finish(counts, sums, reference_codes, allele_sites):
    n = counts.shape[-1]
    features = np.zeros((n, len(CHANNELS)), dtype=np.float32)
    values = {}
    def ratio(a, b):
        return np.divide(a, b, out=np.zeros(n, dtype=np.float64), where=b > 0)
    for group, prefix in enumerate(('', 'hp1_', 'hp2_')):
        c = counts[group]
        total = c[:5].sum(axis=0, dtype=np.uint64)
        if np.any(total > np.iinfo(np.uint32).max):
            raise OverflowError('uint32 depth 범위 초과')
        denominator = total+c[INDEX['deletion']].astype(np.uint64)
        def put(name, array):
            key = prefix+name
            if key in CHANNEL_INDEX:
                features[:, CHANNEL_INDEX[key]] = array
        put('log1p_depth', np.log1p(total.astype(np.float64)))
        if group == 0:
            put('strand_forward_ratio', ratio(c[INDEX['forward']], total))
            put('strand_reverse_ratio', ratio(total-c[INDEX['forward']], total))
        for base, name in enumerate('ACGTN'):
            put('read_'+name+'_ratio', ratio(c[base], total))
        put('insertion_ratio', ratio(c[INDEX['insertion']], denominator))
        put('deletion_after_event_ratio', ratio(c[INDEX['deletion_event']], denominator))
        put('deletion_after_max_len', c[INDEX['max_deletion']])
        values[prefix+'depth'] = total.astype(np.uint32)
        values[prefix+'supplementary_depth'] = c[INDEX['supplementary']].copy()
        values[prefix+'max_deletion_length'] = c[INDEX['max_deletion']].copy()
    total = values['depth'].astype(np.uint64)
    features[:, CHANNEL_INDEX['mean_mapq']] = ratio(sums[0, 0], total)
    features[:, CHANNEL_INDEX['mean_baseq']] = ratio(sums[0, 1], counts[0, INDEX['baseq_n']])
    values['baseq_depth'] = counts[0, INDEX['baseq_n']].copy()
    valid_ref = reference_codes < 4
    positions = np.flatnonzero(valid_ref)
    reference_count = np.zeros(n, dtype=np.uint32)
    reference_count[positions] = counts[0, reference_codes[positions], positions]
    alternatives = counts[0, :4].copy()
    alternatives[reference_codes[positions], positions] = 0
    max_alt = alternatives.max(axis=0)
    features[:, CHANNEL_INDEX['read_supports_variant_ratio']] = np.where(valid_ref, ratio(max_alt, total), 0)
    features[:, CHANNEL_INDEX['base_diff_from_ref_ratio']] = np.where(valid_ref, ratio(total-reference_count, total), 0)
    refs, alts, site, heterozygous = allele_sites
    positions = np.flatnonzero(site)
    ref_count = np.zeros(n, dtype=np.uint32)
    alt_count = np.zeros(n, dtype=np.uint32)
    ref_count[positions] = counts[0, refs[positions], positions]
    alt_count[positions] = counts[0, alts[positions], positions]
    denominator = ref_count.astype(np.uint64)+alt_count
    ref_fraction, baf = ratio(ref_count, denominator), ratio(alt_count, denominator)
    features[:, CHANNEL_INDEX['vcf_ref_allele_fraction']] = ref_fraction
    features[:, CHANNEL_INDEX['baf']] = baf
    features[:, CHANNEL_INDEX['minor_allele_fraction']] = np.minimum(ref_fraction, baf)
    if not np.isfinite(features).all():
        raise ValueError('pileup feature NaN/inf')
    values.update(features=features, reference_base_code=reference_codes,
                  baf_site=site, baf_heterozygous=heterozygous,
                  baf_valid=(site.astype(bool) & (denominator > 0)).astype(np.uint8))
    return values


class ReferenceCache:
    def __init__(self, bam, contig, config, reference, variants, sample, *, oracle=False,
                 block=262144, cache_bytes=512*1024**2):
        if block <= 0 or cache_bytes <= 0:
            raise ValueError('pileup block/cache는 양수여야 함')
        self.bam, self.contig, self.config = bam, contig, config
        self.block, self.limit = block, cache_bytes
        self.length = bam.get_reference_length(contig)
        self.cache = OrderedDict()
        self.bytes = self.blocks_computed = self.cache_hits = 0
        self.native = self.reference = self.variants = None
        self.backend = 'python_cigar_hp'
        try:
            self.reference = pysam.FastaFile(str(reference))
            self.variants = VariantLookup(variants, sample)
            key = 'MFM_PILEUP_VALIDATION_BACKEND' if oracle else 'MFM_PILEUP_PRODUCTION_BACKEND'
            requested = os.environ.get(key, 'auto')
            allowed = {'auto', 'native', 'aligned_pairs'} if oracle else {'auto', 'native', 'python'}
            if requested not in allowed:
                raise ValueError(key+' 설정 오류')
            if requested not in {'python', 'aligned_pairs'}:
                from .complete_native import CompleteNative
                try:
                    self.native = CompleteNative(bam, contig, config, oracle=oracle)
                    self.backend = 'native_htslib_hp' if oracle else 'native_cigar_hp'
                except NativeUnavailable:
                    if requested == 'native' or oracle:
                        raise
            elif oracle:
                raise NativeUnavailable('v2 source 검증은 독립 HTSlib native를 요구함')
        except BaseException:
            self.close()
            raise

    def close(self):
        for handle in (self.native, self.reference, self.variants):
            if handle is not None:
                handle.close()
        self.cache.clear()
        self.bytes = 0

    def _block(self, key):
        value = self.cache.pop(key, None)
        if value is not None:
            self.cache_hits += 1
        else:
            start, end = key*self.block, min((key+1)*self.block, self.length)
            counts, sums = (self.native.calculate(start, end) if self.native is not None else
                            calculate_counts(self.bam, self.contig, start, end, self.config))
            reference = self.reference.fetch(self.contig, start, end).upper().encode('ascii')
            codes = TABLE[np.frombuffer(reference, dtype=np.uint8)].copy()
            codes[codes > 4] = 4
            value = finish(counts, sums, codes, self.variants.block(self.contig, start, end, codes))
            self.blocks_computed += 1
            self.bytes += sum(array.nbytes for array in value.values())
        self.cache[key] = value
        while self.bytes > self.limit and len(self.cache) > 1:
            _, removed = self.cache.popitem(last=False)
            self.bytes -= sum(array.nbytes for array in removed.values())
        return value

    def gather(self, positions, names=None):
        names = names or tuple(name for name in ARRAYS if name != 'mapping_kind')
        out = {name: np.zeros((len(positions), len(CHANNELS)) if name == 'features' else len(positions),
                              dtype=ARRAYS[name]) for name in names}
        if 'reference_base_code' in out:
            out['reference_base_code'].fill(255)
        valid = (positions >= 0) & (positions < self.length)
        for key in np.unique(positions[valid]//self.block):
            value = self._block(int(key))
            selected = valid & (positions//self.block == key)
            index = positions[selected]-int(key)*self.block
            for name in names:
                out[name][selected] = value[name][index]
        return out


def projected(cache, projection, start, actual, padded):
    positions, anchors, kind = projection.window(start, actual, padded)
    result = cache.gather(positions)
    inserted = kind == 2
    if inserted.any():
        anchor = cache.gather(anchors[inserted], ('features',))['features']
        for channel in INSERTION_CHANNELS:
            result['features'][inserted, channel] = anchor[:, channel]
    result['mapping_kind'] = kind
    return result
