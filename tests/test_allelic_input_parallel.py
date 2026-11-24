"""cross-contig fragment·동점·duplicate와 모든 count 출력 동등성을 검증함."""
from pathlib import Path
import importlib.util,json,os,subprocess,sys,tempfile,unittest
from unittest.mock import patch
import pysam
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from integrate_wgs_rna.allelic_input import read_as,read_pair_as


def make_bam(path,index=True):
    header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000},{'SN':'chr2','LN':1000}]}
    with pysam.AlignmentFile(str(path),'wb',header=header) as output:
        for chrom in [0,1]:
            for i,name in enumerate(['cross','tie','dup','multi','no_as','ordinary']):
                r=pysam.AlignedSegment();r.query_name=name;r.query_sequence='A'*20
                r.flag=64 if chrom==0 else 128;r.reference_id=chrom;r.reference_start=10+i*30
                if name=='tie':r.flag=64
                if name=='dup' and chrom==1:r.flag|=1024
                r.mapping_quality=100 if name=='multi' else 255;r.cigar=[(0,20)]
                if name!='no_as':r.set_tag('AS',50 if name=='tie' else 60+chrom)
                output.write(r)
        r=pysam.AlignedSegment();r.query_name='unmapped';r.query_sequence='A'*20;r.flag=4;r.reference_id=-1;r.reference_start=-1;output.write(r)
    if index:pysam.index(str(path))


class ParallelInputTests(unittest.TestCase):
    def test_all_tables_duplicate_and_header_tie_order_match(self):
        with tempfile.TemporaryDirectory() as raw:
            path=Path(raw)/'input.bam';make_bam(path)
            for duplicates in [False,True]:
                expected=[read_as(str(path),255,duplicates) for _ in range(2)]
                actual,plan=read_pair_as(path,path,255,duplicates,workers=3,benchmark=True,scratch=raw)
                self.assertEqual(expected,actual);self.assertTrue(plan['serial_equal'])
                self.assertEqual(actual[0][0]['tie'][0][1],'chr1')
                self.assertEqual(actual[0][0]['cross'][1][1],'chr2')
                self.assertIn('dup',actual[0][3]);self.assertEqual(actual[0][2],1)
                self.assertEqual(list(Path(raw).glob('ase_input_*')),[])

    def test_unindexed_whole_hap_order_and_memory_preflight(self):
        with tempfile.TemporaryDirectory() as raw:
            path=Path(raw)/'input.bam';make_bam(path,index=False)
            actual,plan=read_pair_as(path,path,255,workers=2,scratch=raw)
            self.assertEqual(plan['mode'],'serial_optimized');self.assertEqual(actual[0],read_as(str(path),255))
            pysam.index(str(path))
            with patch.dict(os.environ,{'MFM_RAM_LIMIT_BYTES':'1000000000'}):
                with self.assertRaisesRegex(ValueError,'RAM'):read_pair_as(path,path,255,workers=2,scratch=raw)
            self.assertEqual(list(Path(raw).glob('ase_input_*')),[])


    def test_gc_state_is_restored_on_success_and_failure(self):
        import gc
        with tempfile.TemporaryDirectory() as raw:
            path=Path(raw)/'input.bam';make_bam(path)
            before=gc.isenabled()
            try:
                for enabled in [True,False]:
                    gc.enable() if enabled else gc.disable()
                    read_as(str(path),255)
                    self.assertEqual(gc.isenabled(),enabled)
                    with self.assertRaises((OSError,ValueError)):
                        read_as(str(path)+'.missing',255)
                    self.assertEqual(gc.isenabled(),enabled)
            finally:
                gc.enable() if before else gc.disable()


if __name__=='__main__':unittest.main()
