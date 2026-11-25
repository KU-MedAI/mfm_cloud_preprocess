"""모델 바이트·외부 링크·삭제 대상·trace 최신성의 기존 계약을 확인함."""
from pathlib import Path
import os
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from handoff.cleanup import model_bytes,model_files,package_inventory
from handoff.metrics import _tree_bytes


class PackageInventoryTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.base=Path(self.temp.name);self.root=self.base/'S';self.root.mkdir()

    def tearDown(self):self.temp.cleanup()

    def write(self,name,data):
        path=self.root/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
        return path

    def test_model_selection_matches_independent_files_and_marks_removable_inputs(self):
        files={
            'compact/layout/chr1/pileup.zarr/HAP1/features/0.0':b'a'*11,
            'compact/layout/manifest.json':b'{}',
            'rna_personal/public/HAP1/coverage.zarr/depth/0':b'b'*7,
            'rna_personal/public/allelic_counts.tsv':b'gene\t1\n',
            'rna_personal/public/report.json':b'not a model',
            'rna_personal/public/STAR.bam':b'raw',
            'logs/source.fa':b'raw fasta',
            'metadata/sample.json':b'{}',
        }
        for name,value in files.items():self.write(name,value)
        expected=11+2+7+7
        scan=package_inventory(self.root)
        self.assertEqual(scan['model_bytes'],expected)
        self.assertEqual(model_bytes(self.root),expected)
        self.assertEqual(sum(p.stat().st_size for p in model_files(self.root)),expected)
        self.assertEqual(_tree_bytes(self.root/'rna_personal',rna=True),14)
        self.assertEqual(set(scan['raw_files']),{self.root/'rna_personal/public/STAR.bam',self.root/'logs/source.fa'})
        self.assertEqual(scan['rna_extras'],[self.root/'rna_personal/public/STAR.bam'])
        self.assertEqual(set(scan['zarr_dirs']),{self.root/'compact/layout/chr1/pileup.zarr',self.root/'rna_personal/public/HAP1/coverage.zarr'})

    def test_file_directory_and_broken_links_are_detected_without_following(self):
        self.write('compact/layout/own',b'own')
        outside=self.base/'outside';outside.mkdir();(outside/'large').write_bytes(b'x'*1000)
        links=[self.root/'compact/layout/file',self.root/'compact/layout/directory',self.root/'rna_personal']
        links[0].symlink_to(outside/'large');links[1].symlink_to(outside,target_is_directory=True)
        links[2].symlink_to(self.base/'missing',target_is_directory=True)
        scan=package_inventory(self.root)
        self.assertEqual(set(scan['symlinks']),set(links))
        self.assertEqual(scan['model_bytes'],3)
        self.assertEqual(model_bytes(self.root),3)

    def test_trace_timestamp_remains_layout_specific_and_fresh(self):
        first=self.write('compact/layout/chr1/a',b'a')
        latest=self.write('compact/layout/chr2/b',b'b')
        unrelated=self.write('compact/unrelated',b'x')
        os.utime(first,ns=(1000,1000));os.utime(latest,ns=(2000,2000));os.utime(unrelated,ns=(3000,3000))
        self.assertEqual(package_inventory(self.root)['newest_compact']['layout'],2000)
        os.utime(latest,ns=(4000,4000))
        self.assertEqual(package_inventory(self.root)['newest_compact']['layout'],4000)


if __name__=='__main__':unittest.main()
