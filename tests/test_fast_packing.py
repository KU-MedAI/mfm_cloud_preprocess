"""4bit byte·offset·오류 계약과 bounded batch 경계를 검증함."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'),str(ROOT/'vendor/MFM-DNARNAFM')]
from integrate_wgs_rna.legacy_shared_hap_converter import pack_4bit_acgtnp_chunk, _pack_4bit_acgtnp_chunk_legacy


class PackingTest(unittest.TestCase):
    def same(self, seq, allowed=frozenset('ACGTNP'), policy='to_n'):
        a = _pack_4bit_acgtnp_chunk_legacy(seq, allowed, 23, policy)
        b = pack_4bit_acgtnp_chunk(seq, allowed, 23, policy)
        for x,y in zip(a[:3], b[:3]):
            np.testing.assert_array_equal(x,y)
        self.assertEqual(a[3:], b[3:])

    def test_odd_empty_lowercase_invalid_and_boundaries(self):
        for seq in [[],[''],['A','C','NNN','acgt','P',''],['RX?','t','nP']]:
            self.same(seq)
        self.same(['N','NNN',''],frozenset('N'))
        self.same(['ACGTN'*31]*70000)

    def test_failure_identifies_first_record_and_position(self):
        for seq in [['A','T?'],['R','é'],['A','é']]:
            messages=[]
            for fn in [_pack_4bit_acgtnp_chunk_legacy,pack_4bit_acgtnp_chunk]:
                try:fn(seq,frozenset('ACGTN'),17,'fail')
                except Exception as exc:messages.append((type(exc),str(exc)))
            self.assertEqual(len(messages),2)
            self.assertEqual(messages[0],messages[1])

if __name__ == '__main__':unittest.main()
