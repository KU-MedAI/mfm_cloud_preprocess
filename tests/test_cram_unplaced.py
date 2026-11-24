"""미배치 CRAM 컨테이너가 여러 개일 때 직렬/병렬 QC의 완전성을 검증한다."""
import gzip
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pysam
from integrate_wgs_rna.realign import qc


class CramUnplacedTest(unittest.TestCase):
    def test_multicontainer_cram_parallel_matches_full_scan(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ref=root/'ref.fa';ref.write_text('>chr1\n'+'ACGT'*250+'\n')
            pysam.faidx(str(ref))
            output=root/'reads.cram'
            header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}]}
            with pysam.AlignmentFile(str(output),'wc',header=header,reference_filename=str(ref),
                                     format_options=[b'seqs_per_slice=10',b'slices_per_container=1']) as dest:
                for i in range(25):
                    r=pysam.AlignedSegment();r.query_name=f'mapped{i}';r.flag=0
                    r.reference_id=0;r.reference_start=i;r.mapping_quality=60;r.cigarstring='40M'
                    r.query_sequence='ACGT'*10;r.query_qualities=pysam.qualitystring_to_array('I'*40)
                    r.set_tag('HP',1);dest.write(r)
                for i in range(67):
                    r=pysam.AlignedSegment();r.query_name=f'unplaced{i}';r.flag=4
                    r.reference_id=-1;r.reference_start=-1;r.query_sequence='TGCA'*10
                    r.query_qualities=pysam.qualitystring_to_array('I'*40)
                    r.set_tag('XU',1);dest.write(r)
            pysam.index(str(output))
            index=Path(str(output)+'.crai');original=index.read_bytes()
            with gzip.open(index,'rt') as fh:
                self.assertGreater(sum(row.startswith('-1\t') for row in fh),1)
            config=SimpleNamespace(output_format='cram',sort_threads=1)
            serial=qc.scan_realigned_output(config,output,ref,1,workers=1)
            parallel=qc.scan_realigned_output(config,output,ref,1,workers=2)
            self.assertEqual(serial,parallel)
            self.assertEqual(parallel['primary'],92)
            self.assertEqual(parallel['rescue_primary'],67)
            self.assertEqual(index.read_bytes(),original)


if __name__=='__main__':
    unittest.main()
