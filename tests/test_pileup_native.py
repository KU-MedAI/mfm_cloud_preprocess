"""서로 다른 집계 경로와 BAM/CRAM·고밀도·CIGAR 경계의 수치 동등성을 확인함."""
from pathlib import Path
from dataclasses import replace
import random
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pysam

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src'),str(ROOT/'vendor/MFM-DNARNAFM')]
from integrate_wgs_rna.pileup_features.compute import calculate, ReferenceCache
from integrate_wgs_rna.pileup_features.native_oracle import NativeOracle, NativeProducer
from integrate_wgs_rna.record_stats_parallel import NativeUnavailable
from integrate_wgs_rna.pileup_features.export import filter_config
from test_pileup_feature import bam_file


class NativePileupTest(unittest.TestCase):
    def test_production_fallback_keeps_calculation_errors_visible(self):
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)/'reads.bam';bam_file(path,'ACGT',[('ACGT',[(0,4)],0,60,30)])
            with pysam.AlignmentFile(str(path),'rb') as bam:
                target='integrate_wgs_rna.pileup_features.native_oracle.NativeProducer'
                with patch.dict('os.environ',{'MFM_PILEUP_PRODUCTION_BACKEND':'auto'}):
                    with patch(target,side_effect=NativeUnavailable('시험용 컴파일 불가')):
                        cache=ReferenceCache(bam,'chr1',filter_config('long_read',30))
                        self.assertEqual(cache.backend,'production_cigar')
                        np.testing.assert_array_equal(cache.gather(np.arange(4))['depth'],[1,1,1,1]);cache.close()
                    with patch(target,side_effect=RuntimeError('시험용 입력 오류')):
                        with self.assertRaisesRegex(RuntimeError,'입력 오류'):ReferenceCache(bam,'chr1',filter_config('long_read',30))

    def test_close_releases_reference_cache_arrays(self):
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)/'reads.bam';bam_file(path,'ACGT',[('ACGT',[(0,4)],0,60,30)])
            with pysam.AlignmentFile(str(path),'rb') as bam:
                for oracle in (False,True):
                    cache=ReferenceCache(bam,'chr1',filter_config('long_read',30),oracle=oracle)
                    cache.gather(np.arange(4))
                    self.assertTrue(cache.cache)
                    cache.close()
                    self.assertFalse(cache.cache)
                    if cache.native is not None:self.assertIsNone(cache.native.handle)

    def test_auto_fallback_only_for_unavailable_native(self):
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)/'reads.bam';bam_file(path,'ACGT',[('ACGT',[(0,4)],0,60,30)])
            with pysam.AlignmentFile(str(path),'rb') as bam:
                config=filter_config('long_read',30)
                target='integrate_wgs_rna.pileup_features.native_oracle.NativeOracle'
                with patch(target,side_effect=NativeUnavailable('시험용 빌드 불가')):
                    cache=ReferenceCache(bam,'chr1',config,oracle=True,validation_backend='auto')
                    self.assertEqual(cache.backend,'aligned_pairs')
                    np.testing.assert_array_equal(cache.gather(np.arange(4))['depth'],[1,1,1,1])
                    cache.close()
                    with self.assertRaises(NativeUnavailable):ReferenceCache(bam,'chr1',config,oracle=True,validation_backend='native')
                with patch(target,side_effect=RuntimeError('시험용 입력 손상')):
                    with self.assertRaisesRegex(RuntimeError,'입력 손상'):ReferenceCache(bam,'chr1',config,oracle=True,validation_backend='auto')

    def compare(self,bam,start,end,config):
        native=NativeOracle(bam,'chr1',config)
        try:actual=native.calculate(start,end)
        finally:native.close()
        for oracle in (False,True):
            expected=calculate(bam,'chr1',start,end,config,oracle=oracle)
            for name in expected:
                np.testing.assert_array_equal(actual[name],expected[name],err_msg=f'{oracle}:{name}:{start}-{end}')
        producer=NativeProducer(bam,'chr1',config)
        try:produced=producer.calculate(start,end)
        finally:producer.close()
        for name in actual:np.testing.assert_array_equal(produced[name],actual[name],err_msg='native producer: '+name)
        return actual

    def test_complex_cigar_filters_and_quality(self):
        with tempfile.TemporaryDirectory() as temporary:
            p=Path(temporary)/'reads.bam'
            records=[('ACGT',[(0,4)],0,60,40),('ACTGT',[(0,2),(1,1),(2,2),(0,2)],16,60,30),
                     ('ACATT',[(0,2),(2,2),(1,1),(0,2)],0,60,None),
                     ('TTACGT',[(5,2),(4,2),(0,2),(6,1),(0,2),(5,1)],16,60,20),
                     ('TTACGT',[(1,2),(0,4)],0,60,0),('ACGT',[(2,2),(0,4)],0,60,None),
                     ('ACGT',[(0,2),(3,2),(0,2)],0,60,40),('ACNT',[(0,4)],0,60,255)]
            records += [('ACGT',[(0,4)],flag,mq,30) for flag,mq in [(1024,60),(2048,60),(256,60),(512,60),(0,29),(1,60),(3,60)]]
            bam_file(p,'ACGT'*8,records)
            with pysam.AlignmentFile(str(p),'rb') as bam:
                for mode in ('long_read','short_read'):
                    config=filter_config(mode,30)
                    self.compare(bam,0,32,config)
                    self.compare(bam,0,3,config)
                    self.compare(bam,3,16,config)
                    config=replace(config,require_proper_pair=True)
                    self.compare(bam,0,32,config)

    def test_randomized_valid_cigars_and_block_boundaries(self):
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)/'random.bam';rng=random.Random(42)
            header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':256}]}
            with pysam.AlignmentFile(str(path),'wb',header=header) as out:
                for i in range(100):
                    cigar=[(4,rng.randrange(1,4)),(0,rng.randrange(3,10))]
                    for _ in range(4):cigar.extend([(rng.choice([1,2,3,6]),rng.randrange(1,5)),(rng.choice([0,7,8]),rng.randrange(2,10))])
                    length=sum(n for op,n in cigar if op in (0,1,4,7,8))
                    read=pysam.AlignedSegment();read.query_name='r'+str(i);read.query_sequence=''.join(rng.choices('ACGTN',k=length))
                    read.flag=rng.choice([0,16,1024,2048,256,512]);read.reference_id=0;read.reference_start=0
                    read.mapping_quality=rng.choice([29,30,60]);read.cigartuples=cigar
                    if i%4:read.query_qualities=rng.choices([0,10,20,40,255],k=length)
                    out.write(read)
            pysam.index(str(path))
            with pysam.AlignmentFile(str(path),'rb') as bam:
                for mode in ('long_read','short_read'):
                    config=filter_config(mode,30)
                    for a,b in ((0,17),(17,83),(83,256),(0,256)):
                        self.compare(bam,a,b,config)

    def test_depth_is_not_limited_to_8000(self):
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)/'deep.bam';bam_file(path,'ACGT',[('A',[(0,1)],0,60,40)]*9001)
            with pysam.AlignmentFile(str(path),'rb') as bam:
                values=self.compare(bam,0,4,filter_config('long_read',30))
                self.assertEqual(int(values['depth'][0]),9001)

    def test_cram_and_missing_quality(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);reference=root/'ref.fa';reference.write_text('>chr1\nACGTACGT\n');pysam.faidx(str(reference))
            bam=root/'reads.bam';bam_file(bam,'ACGTACGT',[('ACGT',[(0,4)],0,60,40),('ACTGT',[(0,2),(1,1),(0,2)],16,60,None)])
            cram=root/'reads.cram';pysam.view('-C','-T',str(reference),'-o',str(cram),str(bam),catch_stdout=False);pysam.index(str(cram))
            with pysam.AlignmentFile(str(cram),'rc',reference_filename=str(reference)) as handle:
                self.compare(handle,0,8,filter_config('long_read',30))

    def test_invalid_observed_base_is_not_silently_replaced(self):
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)/'invalid.bam';bam_file(path,'ACGT',[('ACRT',[(0,4)],0,60,30)])
            with pysam.AlignmentFile(str(path),'rb') as bam:
                oracle=NativeOracle(bam,'chr1',filter_config('long_read',30))
                try:
                    with self.assertRaisesRegex(RuntimeError,'unsupported observed base'):oracle.calculate(0,4)
                finally:oracle.close()


if __name__=='__main__':unittest.main()
