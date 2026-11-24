"""병합 결과·QC·활성 파일을 보존하며 임시 조각만 제거하는지 검증함."""
from pathlib import Path
import json,tempfile,time,unittest
import pysam
from integrate_wgs_rna.haplotag_cleanup import cleanup_parallel_shards


class Tests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.bam=self.root/'S.bam';self.tags=self.root/'S.tags.tsv.gz';self.qc=self.root/'S.qc.json'
        header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':100}], 'RG':[{'ID':'S','SM':'S'}]}
        with pysam.AlignmentFile(str(self.bam),'wb',header=header) as output:
            r=pysam.AlignedSegment(output.header);r.query_name='r';r.query_sequence='ACGT';r.query_qualities=pysam.qualitystring_to_array('IIII')
            r.reference_id=0;r.reference_start=1;r.flag=0;r.mapping_quality=60;r.cigarstring='4M';r.set_tag('HP',1);output.write(r)
        pysam.index(str(self.bam));self.tags.write_bytes(b'checked')
        self.folder=self.root/'parallel_haplotag/run';self.folder.mkdir(parents=True)
        self.shard=self.folder/'0000.bam';self.shard.write_bytes(b'duplicate')
        self.log=self.folder/'0000.log';self.log.write_text('근거 로그 보존함')
        self.report=self.folder/'execution.json';self.report.write_text(json.dumps({'sample':'S','status':'success','input_alignment':'source.bam','output_bam':str(self.bam),'finished':time.time()-1}))
        self.payload={'schema':'mfm_haplotag_sample_qc_v1','status':'pass','action':'created_haplotagged_bam','sample':'S','input_alignment':'source.bam','haplotagged_bam':str(self.bam),'downstream_alignment':str(self.bam),'taglist':str(self.tags),'after':{'hp_present':1},'taglist_counts':{'total':1}}
        self.qc.write_text(json.dumps(self.payload))

    def tearDown(self):self.temp.cleanup()

    def test_only_temporary_shards_removed(self):
        before=self.bam.read_bytes();result=cleanup_parallel_shards(self.qc)
        self.assertEqual(result['status'],'complete');self.assertFalse(self.shard.exists())
        self.assertEqual(self.bam.read_bytes(),before);self.assertTrue(self.log.exists());self.assertTrue(self.report.exists())

    def test_open_shard_is_preserved(self):
        with self.shard.open('rb'):
            result=cleanup_parallel_shards(self.qc)
            self.assertEqual(result['status'],'deferred_open_files');self.assertTrue(self.shard.exists())

    def test_failed_qc_preserves_everything(self):
        self.payload['status']='fail';self.qc.write_text(json.dumps(self.payload))
        with self.assertRaises(ValueError):cleanup_parallel_shards(self.qc)
        self.assertTrue(self.shard.exists())

    def test_newer_failed_run_is_preserved(self):
        other=self.folder.parent/'new';other.mkdir();candidate=other/'0000.bam';candidate.write_bytes(b'partial')
        (other/'execution.json').write_text(json.dumps({'sample':'S','status':'failed','input_alignment':'source.bam','finished':time.time()+10}))
        cleanup_parallel_shards(self.qc);self.assertTrue(candidate.exists())

    def test_dry_run_does_not_delete(self):
        result=cleanup_parallel_shards(self.qc,dry_run=True)
        self.assertEqual(len(result['files']),1);self.assertTrue(self.shard.exists())

    def test_missing_index_preserves_everything(self):
        Path(str(self.bam)+'.bai').unlink()
        with self.assertRaises(ValueError):cleanup_parallel_shards(self.qc)
        self.assertTrue(self.shard.exists())


if __name__=='__main__':unittest.main()
