"""통계 cache의 완료계약·무효화·현재 RNA 파일 집계·심볼릭링크 제외를 검증함."""
from pathlib import Path
import json,sys,tempfile,unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from handoff.cleanup import model_bytes
from handoff.metrics import statistics_model_bytes,_tree_bytes,snapshot


class StatisticsCacheTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.layout=self.root/'compact/layout';self.layout.mkdir(parents=True)
        self.manifest=self.layout/'manifest.json';self.write_manifest('complete')
        (self.layout/'chr1').mkdir();(self.layout/'chr1/chunk').write_bytes(b'abc')
        self.rna=self.root/'rna_personal/HAP1';self.rna.mkdir(parents=True)
        (self.rna/'coverage.zarr').mkdir();(self.rna/'coverage.zarr/chunk').write_bytes(b'abcde')
        (self.rna/'reads.bam').write_bytes(b'x'*200)
        (self.root/'rna_personal/counts.tsv').write_bytes(b'counts\n')

    def write_manifest(self,status):
        self.manifest.write_text(json.dumps(dict(status=status,records={'chr1':dict(contig='chr1')})))

    def test_completed_cache_avoids_compact_scan_and_includes_new_rna_tables(self):
        before,hit=statistics_model_bytes(self.root)
        self.assertEqual(before,model_bytes(self.root));self.assertFalse(hit)
        (self.root/'rna_personal/ase_qc.tsv').write_bytes(b'new results\n')
        def counted(path,**kwargs):
            self.assertNotEqual(Path(path),self.root/'compact')
            return _tree_bytes(path,**kwargs)
        with patch('handoff.metrics._tree_bytes',side_effect=counted):
            after,hit=statistics_model_bytes(self.root)
        self.assertTrue(hit);self.assertEqual(after,model_bytes(self.root));self.assertEqual(after-before,12)

    def test_manifest_or_metadata_change_invalidates_cache(self):
        statistics_model_bytes(self.root)
        (self.layout/'chr1/chunk').write_bytes(b'updated')
        (self.layout/'chr1/contig_metadata.json').write_text('{"changed":true}')
        value,hit=statistics_model_bytes(self.root)
        self.assertFalse(hit);self.assertEqual(value,model_bytes(self.root))
        self.write_manifest('running');statistics_model_bytes(self.root)
        (self.layout/'chr1/chunk').write_bytes(b'updated again')
        value,hit=statistics_model_bytes(self.root)
        self.assertFalse(hit);self.assertEqual(value,model_bytes(self.root))

    def test_bad_cache_and_links_match_exact_uncached_counter(self):
        (self.root/'metadata').mkdir();(self.root/'metadata/compact_bytes_cache.json').write_text('invalid')
        (self.layout/'chunk_link').symlink_to(self.layout/'chr1/chunk')
        (self.layout/'directory_link').symlink_to(self.layout/'chr1',target_is_directory=True)
        (self.rna/'counts_link.tsv').symlink_to(self.root/'rna_personal/counts.tsv')
        value,hit=statistics_model_bytes(self.root)
        self.assertFalse(hit);self.assertEqual(value,model_bytes(self.root))

    def test_snapshot_prunes_zarr_and_logs_but_keeps_qc_files(self):
        (self.root/'qc').mkdir();(self.root/'qc/result.json').write_text('{}')
        (self.rna/'coverage.zarr/hidden_qc.json').write_text('{}')
        (self.root/'logs').mkdir();(self.root/'logs/hidden_qc.json').write_text('{}')
        row=snapshot(self.root,'sample','25','RNA')
        self.assertEqual(row['qc_files'],['qc/result.json'])


if __name__=='__main__':unittest.main()
