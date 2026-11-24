import json,os,subprocess,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import pysam
from handoff.reference_iupac import normalized_reference
from handoff import cli,planner
ROOT=Path(__file__).resolve().parents[1]
class ReferenceIupacTests(unittest.TestCase):
 def fixture(self,t,sequence='ACGTRYNNACGT'):
  p=Path(t)/'source.fa';p.write_text('>chr1\n'+sequence+'\n');pysam.faidx(str(p));return p
 def test_only_iupac_changes_and_cache_is_stable(self):
  with tempfile.TemporaryDirectory() as t:
   source=self.fixture(t);original=source.read_bytes();target=Path(t)/'normalized.fa'
   self.assertEqual(normalized_reference(source,target),str(target));before=target.stat().st_mtime_ns
   with pysam.FastaFile(str(target)) as fa:self.assertEqual(fa.fetch('chr1'),'ACGTNNNNACGT');self.assertEqual(fa.get_reference_length('chr1'),12)
   record=json.loads(Path(str(target)+'.source.json').read_text());self.assertEqual(record['changed_bases'],2);self.assertEqual(record['canonical_acgtn_changes'],0)
   normalized_reference(source,target);self.assertEqual(target.stat().st_mtime_ns,before);self.assertEqual(source.read_bytes(),original)
 def test_modified_cache_is_rejected(self):
  with tempfile.TemporaryDirectory() as t:
   source=self.fixture(t);target=Path(t)/'normalized.fa';normalized_reference(source,target);target.write_text('>chr1\nAAAAAAAAAAAA\n')
   with self.assertRaises(ValueError):normalized_reference(source,target)
 def test_unknown_non_iupac_character_is_rejected(self):
  with tempfile.TemporaryDirectory() as t:
   source=self.fixture(t,'ACGTZ');target=Path(t)/'normalized.fa'
   with self.assertRaises(ValueError):normalized_reference(source,target)
   self.assertFalse(target.exists())
 def test_binding_is_explicit_and_default_reference_is_preserved(self):
  row=next(r for r in planner.read_tsv(ROOT/'manifests/리드별/전체_long_read.tsv') if r['sample_id']=='HG03492')
  with tempfile.TemporaryDirectory() as t:
   source=self.fixture(t);refs={'autosomes':str(source),'ncbi':str(source),'broad':str(source)}
   config=dict(data_root=str(Path(t)/'input'),output_root=str(Path(t)/'output'),state_root=str(Path(t)/'state'),run_id='fixture',pipeline_env={'MFM_LONG_READ_OUTPUT_MODE':'sequence_only'})
   with patch('handoff.bundled_assets.panel_directory',return_value=Path(t)):
    _,first=cli.bind_manifest(config,row,[],refs);self.assertEqual(first['vcf_reference_fasta'],str(source))
    config['sample_pipeline_env']={'HG03492':{'MFM_REFERENCE_IUPAC_POLICY':'to_n'}}
    _,second=cli.bind_manifest(config,row,[],refs);self.assertNotEqual(second['vcf_reference_fasta'],str(source));self.assertEqual(second['vcf_reference_fasta'],second['analysis_reference_fasta'])
 def test_real_norm_and_hap_alleles_without_dropping_variant(self):
  with tempfile.TemporaryDirectory() as t:
   source=self.fixture(t);target=Path(t)/'normalized.fa';normalized_reference(source,target)
   vcf=Path(t)/'source.vcf';vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=12>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t5\t.\tN\tC\t.\tPASS\t.\tGT\t1|0\n')
   invalid=subprocess.run(['bcftools','norm','-c','e','-f',str(source),str(vcf)],capture_output=True);self.assertNotEqual(invalid.returncode,0)
   normalized=Path(t)/'valid.vcf.gz';subprocess.run(['bcftools','norm','-c','e','-f',str(target),'-Oz','-o',str(normalized),str(vcf)],check=True,capture_output=True);subprocess.run(['bcftools','index','-t',str(normalized)],check=True)
   for hap,expected in [('1','ACGTCNNNACGT'),('2','ACGTNNNNACGT')]:
    result=subprocess.run(['bcftools','consensus','-f',str(target),'-H',hap,str(normalized)],check=True,capture_output=True,text=True);self.assertEqual(''.join(result.stdout.splitlines()[1:]),expected);self.assertIn('Applied 1 variants' if hap=='1' else 'Applied 0 variants',result.stderr)
if __name__=='__main__':unittest.main()
