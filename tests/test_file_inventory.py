"""최신 변경 검사의 파일 링크 처리를 유지하고 외부 디렉터리는 순회하지 않음."""
from pathlib import Path
import os
import tempfile
import unittest
from tools.file_inventory import file_stats


class InventoryTests(unittest.TestCase):
    def test_same_files_and_mtimes_as_previous_rglob_including_file_links(self):
        with tempfile.TemporaryDirectory() as raw:
            base=Path(raw);root=base/'model';root.mkdir()
            (root/'sub').mkdir();(root/'sub/chunk').write_bytes(b'123')
            outside=base/'outside';outside.mkdir();(outside/'chunk').write_bytes(b'45')
            (root/'file_link').symlink_to(outside/'chunk')
            (root/'dir_link').symlink_to(outside,target_is_directory=True)
            (root/'missing').symlink_to(base/'absent')
            expected={str(p):p.stat().st_mtime_ns for p in root.rglob('*') if p.is_file()}
            actual={str(p):info.st_mtime_ns for p,info in file_stats(root)}
            self.assertEqual(actual,expected)
            self.assertNotIn(str(root/'dir_link/chunk'),actual)
            (root/'sub/new').write_bytes(b'new')
            self.assertIn(str(root/'sub/new'),dict((str(p),s) for p,s in file_stats(root)))

    def test_missing_tree_is_empty(self):
        with tempfile.TemporaryDirectory() as raw:self.assertEqual(list(file_stats(Path(raw)/'missing')),[])
