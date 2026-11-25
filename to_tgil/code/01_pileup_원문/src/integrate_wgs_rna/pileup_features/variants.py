"""개인 sequence VCF에서 biallelic SNV의 allele 정의만 bounded 조회함."""
from pathlib import Path
import numpy as np
import pysam

BASES = {base: i for i, base in enumerate('ACGT')}


class VariantLookup:
    def __init__(self, path, sample):
        self.path = Path(path)
        self.handle = pysam.VariantFile(str(self.path))
        if sample not in self.handle.header.samples:
            self.handle.close()
            raise ValueError('BAF VCF에 요청한 sample이 없음: ' + sample)
        self.sample = sample

    def close(self):
        self.handle.close()

    def block(self, contig, start, end, reference_codes):
        n = end-start
        refs = np.full(n, 255, dtype=np.uint8)
        alts = np.full(n, 255, dtype=np.uint8)
        site = np.zeros(n, dtype=np.uint8)
        het = np.zeros(n, dtype=np.uint8)
        if contig not in self.handle.header.contigs:
            alternate = contig[3:] if contig.startswith('chr') else 'chr'+contig
            if alternate in self.handle.header.contigs:
                raise ValueError('BAF VCF/reference contig 이름 불일치: ' + contig)
            return refs, alts, site, het
        seen = set()
        for record in self.handle.fetch(contig, start, end):
            index = record.start-start
            if not 0 <= index < n:
                continue
            if index in seen:
                refs[index] = alts[index] = 255
                site[index] = het[index] = 0
                continue
            seen.add(index)
            alleles = record.alleles
            gt = record.samples[self.sample].get('GT')
            if (not alleles or len(alleles) != 2 or alleles[0] not in BASES or alleles[1] not in BASES
                    or alleles[0] == alleles[1] or gt is None or len(gt) != 2
                    or any(value not in (0, 1) for value in gt)):
                continue
            ref, alt = BASES[alleles[0]], BASES[alleles[1]]
            if int(reference_codes[index]) != ref:
                raise ValueError(f'BAF VCF REF/reference 서열 불일치: {contig}:{record.pos}')
            refs[index], alts[index], site[index], het[index] = ref, alt, 1, int(gt[0] != gt[1])
        return refs, alts, site, het
