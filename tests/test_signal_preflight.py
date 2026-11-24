import argparse,importlib.util,json,os,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import pysam
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'src')]
name='13_rna_to_fastq.py' if (ROOT/'09_personal_fasta.sh').exists() else '11_rna_to_fastq.py'
spec=importlib.util.spec_from_file_location('rna_fastq_signal_recovery',ROOT/name);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
class RecoveryTests(unittest.TestCase):
 def test_real_rna_preflight_created_without_wgs(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t);bam=p/'rna.bam';header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':100}],'RG':[{'ID':'r','SM':'S1'}]}
   with pysam.AlignmentFile(str(bam),'wb',header=header) as out:
    for start,flag in [(10,99),(30,147)]:
     r=pysam.AlignedSegment();r.query_name='pair';r.query_sequence='A'*20;r.query_qualities=pysam.qualitystring_to_array('I'*20);r.flag=flag;r.reference_id=0;r.reference_start=start;r.mapping_quality=255;r.cigar=((0,20),);r.next_reference_id=0;r.next_reference_start=30 if start==10 else 10;r.set_tag('RG','r');out.write(r)
   pysam.index(str(bam));gtf=p/'genes.gtf';gtf.write_text('chr1\tfixture\texon\t1\t80\t.\t+\t.\tgene_id "G"; transcript_id "T";\n')
   args=argparse.Namespace(sample='S1',rna=str(bam),gtf=gtf,contigs='chr1',samtools='samtools',hap_read_mode='all_to_both',sm_expected='S1',min_paired_fraction=.9,ratio_qc_policy='record_only',sm_placeholders='sm1',require_sm_match=1,preflight_json=str(p/'qc.json'),decode_reference='')
   with patch.dict(os.environ,{'MFM_WGS_PAYLOAD':'sequence_only','MFM_RNA_OUTPUT_MODE':'coverage_only','MFM_CPU_BUDGET':'2','MFM_SAMPLE_CPU_BUDGET':'2'}):
    self.assertTrue(module.prepare_missing_preflight(args))
   d=json.loads((p/'qc.json').read_text());self.assertEqual(d['status'],'pass');self.assertEqual(d['checks']['primary_records'],2);self.assertEqual(d['checks']['paired_fraction'],1.0)
   with patch.dict(os.environ,{'MFM_WGS_PAYLOAD':'full','MFM_RNA_OUTPUT_MODE':'full'}):self.assertFalse(module.prepare_missing_preflight(args))
if __name__=='__main__':unittest.main()
