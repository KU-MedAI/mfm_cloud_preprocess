"""실제 정렬된 BAM/CRAM만 sort를 생략하고 read·flag·태그·순서를 보존함."""
from pathlib import Path
import contextlib
import hashlib
import io
import os
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import pysam
from integrate_wgs_rna import prepare_alignment as prepare


class FastInputTests(unittest.TestCase):
    def fixture(self, root, *, reverse_tie=False, declared=True):
        ref=root/'ref.fa';sequence='ACGT'*256
        ref.write_text('>chr1\n'+sequence+'\n');pysam.faidx(str(ref))
        source=root/'source.bam'
        header={'HD':{'VN':'1.6','SO':'coordinate' if declared else 'unsorted'},
                'SQ':[{'SN':'chr1','LN':len(sequence),'M5':hashlib.md5(sequence.encode()).hexdigest()}],
                'RG':[{'ID':'r','SM':'S'}]}
        records=[(10,0),(10,16),(30,1024),(50,0)]
        if reverse_tie:records[:2]=reversed(records[:2])
        with pysam.AlignmentFile(source,'wb',header=header) as out:
            for i,(position,flag) in enumerate(records):
                read=pysam.AlignedSegment();read.query_name='read'+str(i)
                read.query_sequence=sequence[position:position+40];read.query_qualities=[30]*40
                read.flag=flag;read.reference_id=0;read.reference_start=position
                read.mapping_quality=60;read.cigarstring='40M'
                read.set_tags([('RG','r'),('HP',i%2+1),('PS',7)]);out.write(read)
        return source,ref

    def run_case(self, *, cram=False, reverse_tie=False, declared=True):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);source,reference=self.fixture(root,reverse_tie=reverse_tie,declared=declared)
            if cram:
                converted=root/'source.cram'
                pysam.view('-C','-T',str(reference),'-o',str(converted),str(source),catch_stdout=False)
                source=converted
            baseline=root/'baseline.bam'
            subprocess.run([shutil.which('samtools'),'sort','--reference',str(reference),'-o',str(baseline),str(source)],check=True,capture_output=True)
            output=root/'prepared.bam';messages=io.StringIO()
            with (root/'prepare.log').open('wb') as log, \
                 patch.dict(os.environ,{'MFM_RECORD_STATS_BACKEND':'native','MFM_CPU_BUDGET':'4',
                                        'MFM_CPU_POOL_GRANTED':'4','MFM_SAMPLE_CPU_BUDGET':'4'}), \
                 contextlib.redirect_stdout(messages):
                prepare.transcode_aligned_source(dict(path=str(source),format='cram' if cram else 'bam'),
                    output,samtools=shutil.which('samtools'),sort_threads=1,sort_memory='64M',
                    fallback_reference=reference,log=log)
            with pysam.AlignmentFile(output,'rb') as actual,pysam.AlignmentFile(baseline,'rb') as expected:
                self.assertEqual([r.to_string() for r in actual],[r.to_string() for r in expected])
            self.assertEqual('"sort_skipped": true' in messages.getvalue(),not reverse_tie and declared)

    def test_sorted_bam_and_cram_keep_exact_records(self):
        for cram in (False,True):
            with self.subTest(cram=cram):self.run_case(cram=cram)

    def test_coordinate_header_cannot_hide_reverse_tie_disorder(self):
        self.run_case(reverse_tie=True)

    def test_undeclared_order_keeps_existing_sort(self):
        self.run_case(declared=False)

    def test_source_change_after_order_check_is_rejected_before_writing(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);source,reference=self.fixture(root);original=prepare.tag_record_stats
            def changed(*args,**kwargs):
                counts=original(*args,**kwargs)
                info=source.stat();os.utime(source,ns=(info.st_atime_ns,info.st_mtime_ns+1))
                return counts
            with (root/'log').open('wb') as log, \
                 patch.dict(os.environ,{'MFM_RECORD_STATS_BACKEND':'native','MFM_CPU_BUDGET':'4',
                                        'MFM_CPU_POOL_GRANTED':'4','MFM_SAMPLE_CPU_BUDGET':'4'}), \
                 patch.object(prepare,'tag_record_stats',side_effect=changed):
                with self.assertRaisesRegex(prepare.PrepareAlignmentError,'입력이 변경'):
                    prepare.transcode_aligned_source(dict(path=str(source),format='bam'),root/'output.bam',
                        samtools=shutil.which('samtools'),sort_threads=1,sort_memory='64M',
                        fallback_reference=reference,log=log)
            self.assertFalse((root/'output.bam').exists())
