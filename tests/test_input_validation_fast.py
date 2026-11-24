"""전수 카운트·누락 검출·REF cache 무효화를 검사함."""
import importlib.util,json,os,shutil,subprocess,tempfile,unittest
from pathlib import Path
import pysam
from unittest.mock import patch
from integrate_wgs_rna.input_parallel import cpu_budget
from integrate_wgs_rna.alignment_record_stats import native_record_stats
from integrate_wgs_rna.input_preflight import reference_m5,validate_vcfs
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('prepare_reference',ROOT/'src/integrate_wgs_rna/prepare_alignment.py')
prepare=importlib.util.module_from_spec(spec);spec.loader.exec_module(prepare)
SAMTOOLS=shutil.which('samtools');BCFTOOLS=shutil.which('bcftools')


class InputValidationFast(unittest.TestCase):
    def test_actual_grant_limits_requested_cpu(self):
        with patch.dict(os.environ,{'MFM_CPU_BUDGET':'72','MFM_CPU_POOL_GRANTED':'3'}):
            self.assertEqual(cpu_budget(),3)

    def test_parallel_budget_retains_worker_memory_and_releases_after_join(self):
        from integrate_wgs_rna.input_parallel import parallel_calls
        from integrate_wgs_rna.memory_policy import GIB
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'leases').mkdir();(root/'waiters').mkdir();lease=root/'leases/own.json'
            lease.write_text(json.dumps({'cpus':4,'total_cpus':4,'memory_bytes':64*GIB}))
            def task(n):
                row=json.loads(lease.read_text())
                self.assertGreaterEqual(row['memory_bytes'],8*GIB)
                return n*n
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_RAM_LIMIT_BYTES':'1000000000000'}):
                self.assertEqual(parallel_calls(task,[(i,) for i in range(8)],4,memory_per_task=GIB),[i*i for i in range(8)])
            row=json.loads(lease.read_text());self.assertEqual(row['cpus'],1);self.assertEqual(row['memory_bytes'],4*GIB)

    def test_reused_bam_preserves_and_rechecks_source_counters(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ref=root/'ref.fa';ref.write_text('>chr1\n'+'A'*1000+'\n');pysam.faidx(str(ref));source=root/'source.bam'
            with pysam.AlignmentFile(source,'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}],'RG':[{'ID':'r','SM':'S'}]}) as out:
                r=pysam.AlignedSegment();r.query_name='r';r.flag=0;r.reference_id=0;r.reference_start=0;r.query_sequence='AA';r.query_qualities=pysam.qualitystring_to_array('IH');r.cigarstring='2M';r.set_tag('RG','r');r.set_tag('HP',1);out.write(r)
            contract={'sample_id':'S','platform':'ont','input_mode':'legacy_alignment','source':{'alignment':str(source),'files':[],'alignment_states':['aligned']},'references':{'analysis':str(ref),'source_decode':str(ref)}}
            qc=root/'qc.json';kwargs=dict(output=root/'out.bam',qc_json=qc,minimap2=shutil.which('minimap2'),samtools=SAMTOOLS,align_threads=1,sort_threads=1,sort_memory='16M',extract_threads=1,overwrite=False,keep_run_bams=False)
            first=prepare.prepare(contract,**kwargs);second=prepare.prepare(contract,**kwargs)
            self.assertEqual(second['status'],'reused_canonical_bam')
            self.assertEqual(first['source_record_tag_counts'],second['source_record_tag_counts'])
            self.assertEqual(second['source_record_tag_counts'][0],second['output_record_tag_counts'])
            changed=json.loads(qc.read_text());changed['source_record_tag_counts'][0]['mapped_records']=0;qc.write_text(json.dumps(changed))
            with self.assertRaises(prepare.PrepareAlignmentError):prepare.prepare(contract,**kwargs)

    def test_all_flag_tag_and_unmapped_counts_match_including_zero_tags(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);bam=root/'mixed.bam';ref=root/'ref.fa';ref.write_text('>chr1\n'+'A'*10000+'\n');pysam.faidx(str(ref))
            with pysam.AlignmentFile(bam,'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':10000}]}) as out:
                for i,flag in enumerate([0,4,256,2048,2304,1024,1536,3072,3840]):
                    r=pysam.AlignedSegment();r.query_name=str(i);r.query_sequence='A'*10;r.query_qualities=pysam.qualitystring_to_array('I'*10);r.flag=flag;r.reference_id=0;r.reference_start=i*10;r.cigarstring='10M'
                    if i%2:r.set_tag('HP',0)
                    if i%3:r.set_tag('PS',0)
                    out.write(r)
                r=pysam.AlignedSegment();r.query_name='tail';r.flag=4;r.query_sequence='C';r.query_qualities=pysam.qualitystring_to_array('I');r.set_tag('HP',1);out.write(r)
            expected=prepare._tag_record_stats_python(bam)
            self.assertEqual(native_record_stats(bam,samtools=SAMTOOLS),expected)
            cram=root/'mixed.cram';subprocess.run([SAMTOOLS,'view','-C','-T',str(ref),'-o',str(cram),str(bam)],check=True,stderr=subprocess.PIPE)
            self.assertEqual(native_record_stats(cram,ref,samtools=SAMTOOLS),expected)
            # EOF만 정상이어도 내부 손상은 전수 해독에서 실패해야 함.
            content=bam.read_bytes();broken=root/'bad.bam';broken.write_bytes(content[:40]+b'\0'*50+content[90:])
            with self.assertRaises(RuntimeError):native_record_stats(broken,samtools=SAMTOOLS)

    def test_reference_m5_preserves_case_and_cache(self):
        import hashlib
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ref=root/'ref.fa';ref.write_text('>chr1\nAaCNN\n>chr2\nTGCA\n');pysam.faidx(str(ref));out=root/'result.tsv'
            self.assertFalse(reference_m5(ref,['chr1','chr2'],out,root/'cache',2))
            self.assertIn(hashlib.md5(b'AaCNN').hexdigest(),out.read_text())
            self.assertTrue(reference_m5(ref,['chr1','chr2'],out,root/'cache',2))
            ref.write_text('>chr1\nAaCNA\n>chr2\nTGCA\n')
            self.assertFalse(reference_m5(ref,['chr1','chr2'],out,root/'cache',2))
            self.assertIn(hashlib.md5(b'AaCNA').hexdigest(),out.read_text())

    def test_ref_mismatch_and_changed_input_do_not_reuse_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ref=root/'ref.fa';ref.write_text('>chr1\n'+'A'*100+'\n');pysam.faidx(str(ref));vcf=root/'x.vcf';gz=root/'x.vcf.gz'
            def write(allele):
                vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t2\t.\t'+allele+'\tT\t.\tPASS\t.\tGT\t0|1\t1|1\n')
                subprocess.run([BCFTOOLS,'view','-Oz','-o',str(gz),str(vcf)],check=True)
                subprocess.run([BCFTOOLS,'index','-f','-t',str(gz)],check=True)
            write('A');args=([(str(gz),'panel_chr1')],str(ref),'chr1',BCFTOOLS,root/'cache',root/'errors','S',1)
            first=validate_vcfs(*args);self.assertTrue(first[0][1]);self.assertFalse(first[0][3])
            second=validate_vcfs(*args);self.assertTrue(second[0][1]);self.assertTrue(second[0][3])
            write('C');third=validate_vcfs(*args);self.assertFalse(third[0][1]);self.assertFalse(third[0][3])
            self.assertIn('Reference allele mismatch',third[0][2])
            Path(str(gz)+'.tbi').unlink();self.assertFalse(validate_vcfs(*args)[0][1])


if __name__=='__main__':unittest.main()
