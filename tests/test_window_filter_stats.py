"""실제 창 생성·구형 캐시·최종 입력 연결에서 제외 사유와 유지 창을 검증한다."""
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

from integrate_wgs_rna.window_index import WindowIndexConfig, run_window_index_generation
from pre_training.preprocess_raw_readset.shared_hap_writer import _apply_external_window_index
from pre_training.preprocess_raw_readset.window_filter_stats import count_index_exclusions
from pre_training.preprocess_raw_readset.windowing import make_windows


class WindowFilterStatsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.seq = list('ACGT' * 100)
        self.seq[150] = self.seq[250] = 'N'
        self.seq = ''.join(self.seq)
        for hap in ('HAP1', 'HAP2'):
            path = self.root/'S1/individual_fasta'/hap/'chr1.fa'
            path.parent.mkdir(parents=True)
            path.write_text('>chr1\n'+self.seq+'\n')
            pysam.faidx(str(path))
        self.bed = self.root/'blacklist.bed'
        self.bed.write_text('chr1\t0\t20\nchr1\t100\t120\nchr1\t300\t301\n')
        self.cfg = WindowIndexConfig(cohort='TEST',output_root=self.root,samples=['S1'],chroms=['chr1'],
                                     sequence_variant_policy='snp_only',window_size=100,offset_unit=1,
                                     blacklist_bed=self.bed,jobs=1,drop_n=True,max_blacklist_frac=.01)
        self.tsv = self.root/'S1/windows/windows.tsv'
        self.meta = self.tsv.with_name('windows.json')

    def read_rows(self):
        with self.tsv.open() as fh:
            return list(csv.DictReader(fh,delimiter='\t'))

    def write_rows(self, rows):
        with self.tsv.open('w',newline='') as fh:
            w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
            w.writeheader();w.writerows(rows)

    def final_selection(self):
        windows=make_windows(self.seq,source_start=0,offset=0,window_size=100,max_n_ratio=0)
        return _apply_external_window_index(window_index_path=self.tsv,label='HAP1',contig='chr1',
                                            coverage_rows=[(w.actual_length,1.,w) for w in windows])

    def test_window_generation_and_final_selection_agree(self):
        result=run_window_index_generation(self.cfg)
        self.assertEqual(result['fail_count'],0)
        meta=json.loads(self.meta.read_text())
        self.assertEqual(meta['drop_reason_priority'],'blacklist_then_n')
        self.assertEqual(meta['n_fraction_denominator'],'total_windows_before_filters')
        for hap in ('HAP1','HAP2'):
            stats=meta['haplotypes'][hap]
            self.assertEqual((stats['total_windows'],stats['dropped_blacklist'],stats['dropped_n'],stats['kept_windows']),(4,2,1,1))
        rows=[r for r in self.read_rows() if r['haplotype']=='HAP1']
        self.assertEqual([r['drop_reason'] for r in rows],['blacklist','blacklist','n',''])
        selected,summary=self.final_selection()
        self.assertEqual([r[2].raw_start for r in selected],[300])
        self.assertEqual(summary['window_index_dropped_blacklist'],2)
        self.assertEqual(summary['window_index_dropped_n'],1)
        self.assertEqual(summary['window_index_drop_reason_priority'],'blacklist_then_n')
        # 마지막 창의 blacklist 1bp/100bp는 정확히1%이므로 유지해야 한다.
        self.assertEqual(rows[-1]['bl_frac'],'0.010000')

    def test_legacy_cache_updates_counts_without_recreating_windows(self):
        run_window_index_generation(self.cfg)
        rows=self.read_rows()
        for r in rows:
            del r['drop_reason']
        self.write_rows(rows)
        meta=json.loads(self.meta.read_text())
        del meta['drop_reason_priority']
        for stats in meta['haplotypes'].values():
            stats['dropped_n'],stats['dropped_blacklist']=2,1
        self.meta.write_text(json.dumps(meta))
        original=self.tsv.read_bytes()
        for path in (self.root/'S1/individual_fasta').rglob('*.fa*'):
            path.unlink()
        with patch('integrate_wgs_rna.window_index.windows_for_fasta',side_effect=AssertionError('FASTA 재계산 금지')):
            result=run_window_index_generation(self.cfg)
        self.assertEqual(result['skip_count'],1)
        self.assertEqual(self.tsv.read_bytes(),original)
        for stats in json.loads(self.meta.read_text())['haplotypes'].values():
            self.assertEqual((stats['dropped_blacklist'],stats['dropped_n']),(2,1))
        selected,summary=self.final_selection()
        self.assertEqual([r[2].raw_start for r in selected],[300])
        self.assertEqual((summary['window_index_dropped_blacklist'],summary['window_index_dropped_n']),(2,1))

    def test_legacy_final_counts_use_original_custom_threshold(self):
        self.tsv.parent.mkdir(parents=True,exist_ok=True)
        self.write_rows([
            dict(contig='chr1',haplotype='HAP1',raw_start=0,actual_length=100,n_nonacgt=1,bl_frac='0.200000',kept=0),
            dict(contig='chr1',haplotype='HAP1',raw_start=100,actual_length=100,n_nonacgt=1,bl_frac='0.400000',kept=0),
            dict(contig='chr1',haplotype='HAP1',raw_start=300,actual_length=100,n_nonacgt=0,bl_frac='0.000000',kept=1),
        ])
        self.meta.write_text(json.dumps(dict(max_blacklist_frac=.3,drop_n=True)))
        _,summary=self.final_selection()
        self.assertEqual((summary['window_index_dropped_blacklist'],summary['window_index_dropped_n']),(1,1))

    def test_exact_decision_survives_rounded_blacklist_fraction(self):
        row=dict(n_nonacgt='1',bl_frac='0.010000',kept='0')
        with self.assertRaisesRegex(ValueError,'반올림 경계'):
            count_index_exclusions([row],max_blacklist_frac=.01)
        for reason,counts in [('blacklist',dict(dropped_blacklist=1,dropped_n=0)),('n',dict(dropped_blacklist=0,dropped_n=1))]:
            self.assertEqual(count_index_exclusions([dict(row,drop_reason=reason)],max_blacklist_frac=.01),counts)
        self.assertEqual(count_index_exclusions([row],max_blacklist_frac=.01,drop_n=False),dict(dropped_blacklist=1,dropped_n=0))


if __name__=='__main__':
    unittest.main()
