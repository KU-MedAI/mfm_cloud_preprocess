"""실제 Stage01의 RNA 보류와 필수 입력 검증을 검사함."""
from pathlib import Path
import os
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]


class Stage01RnaDeferTests(unittest.TestCase):
    def test_deferred_modes_skip_missing_rna_but_required_modes_still_fail(self):
        with tempfile.TemporaryDirectory() as raw:
            base=Path(raw);missing=base/'not_downloaded_rna.bam'
            env=dict(os.environ,PIPELINE_MANIFEST='',RNA_PATH=str(missing),
                DEFAULT_PLATFORM='hifi',DEFAULT_READ_TYPE='long_read',OUTPUT_ROOT=str(base/'output'),
                OUTPUT_ROOT_OVERRIDE='1',INPUT_ROOT=str(base/'input'),
                PYTHON=sys.executable,PYTHON_BIN=sys.executable,
                MFM_CORE_ENV=str(Path(sys.executable).parent.parent),
                MFM_CPU_BUDGET='2',MFM_SAMPLE_CPU_BUDGET='2',
                MFM_WGS_BACKEND='cpu',MFM_RNA_BACKEND='cpu',MFM_RNA_OUTPUT_MODE='full')
            for mode,deferred,success in [('pileup_feature','1',True),('sequence_only','1',True),
                                          ('pileup_feature','0',False),('full','1',False)]:
                with self.subTest(mode=mode,deferred=deferred):
                    result=subprocess.run(['bash',str(ROOT/'01_prepare_analysis_bam.sh'),'S1'],cwd=ROOT,
                        env=dict(env,MFM_LONG_READ_OUTPUT_MODE=mode,MFM_SEQUENCE_DEFER_RNA=deferred),
                        capture_output=True,text=True,timeout=20)
                    self.assertEqual(result.returncode==0,success,result.stdout+result.stderr)
                    if success:self.assertIn('legacy 실행',result.stdout+result.stderr)
                    else:self.assertIn(str(missing),result.stdout+result.stderr)
            self.assertFalse(missing.exists())

    def test_stage02_keeps_wgs_contract_and_defers_only_optional_rna(self):
        import pysam
        with tempfile.TemporaryDirectory() as raw:
            base=Path(raw);reference=base/'ref.fa';bam=base/'reads.bam';vcf=base/'variants.vcf'
            reference.write_text(''.join(f'>chr{i}\n'+('A'*100)+'\n' for i in range(1,23)))
            pysam.faidx(str(reference))
            header={'HD':{'VN':'1.6','SO':'coordinate'},
                    'SQ':[{'SN':f'chr{i}','LN':100,'M5':hashlib.md5(b'A'*100).hexdigest()} for i in range(1,23)],
                    'RG':[{'ID':'rg1','SM':'S1'}]}
            with pysam.AlignmentFile(str(bam),'wb',header=header) as out:
                read=pysam.AlignedSegment();read.query_name='read1';read.query_sequence='A'*10
                read.flag=0;read.reference_id=0;read.reference_start=0;read.mapping_quality=60
                read.cigarstring='10M';read.query_qualities=pysam.qualitystring_to_array('I'*10)
                read.set_tag('RG','rg1');out.write(read)
            pysam.index(str(bam))
            vcf.write_text('##fileformat=VCFv4.2\n'+''.join(f'##contig=<ID=chr{i},length=100>\n' for i in range(1,23))+
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'+
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'+
                'chr1\t2\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\n')
            pysam.tabix_compress(str(vcf),str(vcf)+'.gz');pysam.tabix_index(str(vcf)+'.gz',preset='vcf')
            env=dict(os.environ,PIPELINE_MANIFEST='',RNA_PATH=str(base/'missing_rna.bam'),
                DEFAULT_PLATFORM='hifi',DEFAULT_READ_TYPE='long_read',DEFAULT_DATASET_ROLE='pre_training',
                DEFAULT_PROJECT_ID='fixture',DEFAULT_COHORT='fixture',ALIGNMENT_PATH=str(bam),
                PREPARED_ALIGNMENT_PATH=str(bam),VCF_SOURCE_POLICY='manifest',
                LONGREAD_SEQUENCE_VCF_SOURCE='public',LONGREAD_HAPLOTAG_VCF_SOURCE='public',
                PUBLIC_VCF_PATH=str(vcf)+'.gz',PUBLIC_VCF_SAMPLE_ID='S1',VCF_REFERENCE_FASTA=str(reference),
                VCF_CALLSET='fixture',PHASE_SOURCE='provided',REF_DECODE=str(reference),REF_SEQUENCE=str(reference),
                PYTHON=sys.executable,PYTHON_BIN=sys.executable,MFM_CORE_ENV=str(Path(sys.executable).parent.parent),
                MFM_CPU_BUDGET='2',MFM_SAMPLE_CPU_BUDGET='2',MFM_WGS_BACKEND='cpu',MFM_RNA_BACKEND='cpu',
                MFM_RNA_OUTPUT_MODE='full',RNA_INDEX_ADDNL_SEQNS='0',PREFLIGHT_MIN_FREE_GIB='0')
            for mode,deferred,success in [('pileup_feature','1',True),('sequence_only','1',True),
                                          ('pileup_feature','0',False),('full','1',False)]:
                with self.subTest(mode=mode,deferred=deferred):
                    run=f'{mode}_{deferred}';output=base/run
                    result=subprocess.run(['bash',str(ROOT/'02_preflight.sh'),'S1'],cwd=ROOT,
                        env=dict(env,DATASET_ROOT=str(output),OUTPUT_ROOT=str(output),OUTPUT_ROOT_OVERRIDE='1',
                            INPUT_ROOT=str(base/'input'),RUN_ID=run,MFM_LONG_READ_OUTPUT_MODE=mode,
                            MFM_SEQUENCE_DEFER_RNA=deferred),capture_output=True,text=True,timeout=40)
                    self.assertEqual(result.returncode==0,success,result.stdout+result.stderr)
                    report=list(output.rglob('S1.preflight.json'));self.assertEqual(len(report),1)
                    checks={row['check']:row for row in json.loads(report[0].read_text())['checks']}
                    self.assertEqual(checks['public_vcf_ref_alleles']['status'],'PASS')
                    self.assertEqual(checks['rna_input_contract']['status'],'PASS' if success else 'FAIL')
                    if success:self.assertEqual(checks['rna_input_contract']['detail'],'DEFERRED')


if __name__=='__main__':unittest.main()
