"""빠른 중간 CRAM이 정렬·품질·태그를 동일하게 복원하는지 확인함."""
from pathlib import Path
from types import SimpleNamespace
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pysam

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src'),str(ROOT/'vendor/MFM-DNARNAFM')]
from integrate_wgs_rna.personal_realign import cram_output_options


class CramOutputTests(unittest.TestCase):
    def test_default_explicit_and_invalid_levels(self):
        config=SimpleNamespace(output_format='cram')
        with patch.dict(os.environ,{},clear=True):
            self.assertIn('cram,level=1',cram_output_options(config,Path('ref with space.fa')))
            self.assertEqual(shlex.split(cram_output_options(config,Path('ref with space.fa')))[-1],'ref with space.fa')
        with patch.dict(os.environ,{'MFM_WGS_CRAM_LEVEL':'5'}):
            self.assertIn('cram,level=5',cram_output_options(config,Path('ref.fa')))
        for bad in ['-1','10','auto','1; touch bad']:
            with patch.dict(os.environ,{'MFM_WGS_CRAM_LEVEL':bad}):
                with self.assertRaisesRegex(ValueError,'0–9'):cram_output_options(config,Path('ref.fa'))

    def test_lossless_records_at_both_compression_levels(self):
        samtools=shutil.which('samtools')
        self.assertIsNotNone(samtools,'이미지/환경에 samtools가 필요함')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ref=root/'ref.fa';ref.write_text('>chr1\n'+'ACGT'*32+'\n');pysam.faidx(str(ref))
            bam=root/'source.bam'
            with pysam.AlignmentFile(str(bam),'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':128}]}) as output:
                for i in range(3):
                    read=pysam.AlignedSegment();read.query_name='read'+str(i);read.query_sequence='ACGT'*5
                    read.reference_id=0;read.reference_start=4+i*4;read.cigarstring='20M';read.mapping_quality=60-i
                    read.flag=[0,16,1024][i]
                    if i<2:read.query_qualities=list(range(20,40))
                    read.set_tag('HP',i%2+1);read.set_tag('PS',123);read.set_tag('ZX','original')
                    output.write(read)
            records=[]
            for level in ['5','1']:
                path=root/(level+'.cram')
                with patch.dict(os.environ,{'MFM_WGS_CRAM_LEVEL':level}):
                    options=shlex.split(cram_output_options(SimpleNamespace(output_format='cram'),ref))
                subprocess.run([samtools,'sort',*options,'-o',str(path),str(bam)],check=True,capture_output=True)
                subprocess.run([samtools,'index',str(path)],check=True,capture_output=True)
                with pysam.AlignmentFile(str(path),'rc',reference_filename=str(ref)) as handle:
                    records.append([(r.query_name,r.flag,r.reference_start,r.cigarstring,r.mapping_quality,r.query_sequence,
                                     list(r.query_qualities) if r.query_qualities is not None else None,sorted(r.get_tags()))
                                    for r in handle.fetch(until_eof=True)])
            self.assertEqual(records[0],records[1])
            self.assertEqual([r[0] for r in records[1]],['read0','read1','read2'])


if __name__=='__main__':unittest.main()
