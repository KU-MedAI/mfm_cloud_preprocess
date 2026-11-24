#!/usr/bin/env python3
"""ASE fragment 계수에서 현행·구번호가 공유하는 좌표/분할 helper.

Stage 21은 validated-HAP 증거 정책이 추가되어 할당 본체를 별도로 유지한다.
Stage 16/20의 기존 할당 함수는 역사적 결과 재현을 위해 legacy 이름으로 보존한다.
"""
from __future__ import annotations

import collections
import hashlib

MATES = (0, 1)  # 0 = read1, 1 = read2 또는 unpaired
EXON_BIN = 16384


def split_bucket(query_name: str) -> int:
    """재실행·프로세스에 무관한 결정적 fragment split을 만든다."""
    digest = hashlib.blake2b(query_name.encode("utf-8"), digest_size=8,
                             person=b"mfm-ase").digest()
    return digest[0] & 1


def load_exons(path):
    ex = collections.defaultdict(list)
    for line in open(path):
        gene, contig, start, end = line.rstrip("\n").split("\t")
        ex[contig].append((int(start), int(end), gene))
    for contig in ex:
        ex[contig].sort()
    return ex


def gene_set(ex):
    return {gene for items in ex.values() for _, _, gene in items}


def restrict_to_common(ex1, ex2):
    """구번호 Stage 16/20의 동일한 양쪽 HAP 유전자 교집합 처리."""
    common = gene_set(ex1) & gene_set(ex2)
    dropped = len(gene_set(ex1) | gene_set(ex2)) - len(common)
    out = []
    for ex in (ex1, ex2):
        filt = collections.defaultdict(list)
        for contig, items in ex.items():
            keep = [(start, end, gene) for start, end, gene in items if gene in common]
            if keep:
                filt[contig] = keep
        out.append(filt)
    return out[0], out[1], dropped


def genes_for(ex, contig, blocks):
    """선형 검색 기준 구현. 인덱스 구현과의 결과 대조용이다."""
    hits = set()
    items = ex.get(contig)
    if not items:
        return hits
    for block_start, block_end in blocks:
        for start, end, gene in items:
            if start >= block_end:
                break
            if end > block_start:
                hits.add(gene)
    return hits


def build_exon_index(ex, bin_size: int = EXON_BIN):
    """(contig, fixed-width bin)에서 겹치는 exon을 빠르게 찾도록 색인한다."""
    index: dict[tuple[str, int], list] = {}
    for contig, items in ex.items():
        for start, end, gene in items:
            last = (end - 1) // bin_size if end > start else start // bin_size
            for bin_id in range(start // bin_size, last + 1):
                index.setdefault((contig, bin_id), []).append((start, end, gene))
    return index


def genes_for_indexed(idx, contig, blocks, bin_size: int = EXON_BIN):
    hits = set()
    for block_start, block_end in blocks:
        if block_end <= block_start:
            continue
        for bin_id in range(block_start // bin_size, (block_end - 1) // bin_size + 1):
            for start, end, gene in idx.get((contig, bin_id), ()):
                if start < block_end and end > block_start:
                    hits.add(gene)
    return hits


def genes_of(idx, slot):
    """해당 fragment의 두 mate가 겹치는 유전자 합집합."""
    hits = set()
    for mate_id in MATES:
        mate = slot[mate_id]
        if mate is not None:
            hits |= genes_for_indexed(idx, mate[1], mate[2])
    return hits


def assign_fragments(A, B, ex1, ex2, require_both_haps: bool):
    """구번호 Stage 16/20의 AS 비교·분할 카운트 규칙(현행 Stage 21은 별도)."""
    a_counts = collections.Counter()
    b_counts = collections.Counter()
    both_a = collections.Counter()
    both_b = collections.Counter()
    split_a = (collections.Counter(), collections.Counter())
    split_b = (collections.Counter(), collections.Counter())
    stats = dict(tie=0, only1=0, only2=0, assigned=0, nocommon=0, dropped_single=0)
    for key in set(A) | set(B):
        ra, rb = A.get(key), B.get(key)
        if ra is not None and rb is not None:
            common = [mate_id for mate_id in MATES if ra[mate_id] is not None and rb[mate_id] is not None]
            if not common:
                stats["nocommon"] += 1
                continue
            score_a = sum(ra[mate_id][0] for mate_id in common)
            score_b = sum(rb[mate_id][0] for mate_id in common)
            if score_a == score_b:
                stats["tie"] += 1
                continue
            bucket = split_bucket(key)
            if score_a > score_b:
                for gene in genes_of(ex1, ra):
                    a_counts[gene] += 1
                    split_a[bucket][gene] += 1
                    both_a[gene] += 1
            else:
                for gene in genes_of(ex2, rb):
                    b_counts[gene] += 1
                    split_b[bucket][gene] += 1
                    both_b[gene] += 1
        else:
            if require_both_haps:
                stats["dropped_single"] += 1
                continue
            bucket = split_bucket(key)
            if ra is not None:
                stats["only1"] += 1
                for gene in genes_of(ex1, ra):
                    a_counts[gene] += 1
                    split_a[bucket][gene] += 1
            else:
                stats["only2"] += 1
                for gene in genes_of(ex2, rb):
                    b_counts[gene] += 1
                    split_b[bucket][gene] += 1
        stats["assigned"] += 1
    return dict(a=a_counts, b=b_counts, both_a=both_a, both_b=both_b,
                split_a=split_a, split_b=split_b), stats
