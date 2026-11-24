"""외부 자산 폴더 없이 공통 입력을 복원하고 재사용하는지 검사한다."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import redirect_stdout
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pysam
from handoff import asset_restore, bundled_assets, download_only, planner


class AssetAutoPrepareTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.base = self.root / 'assets'
        self.config = {'bundled_assets_root': str(self.base),
                       'data_root': str(self.root / 'input')}
        self.source = self.root / 'source.bin'
        self.source.write_bytes(b'fixture-panel')
        self.relative = 'common_vcf/panel.vcf.gz'
        self.catalog = self.root / 'manifests/공통자산.json'
        self.catalog.parent.mkdir()
        self.write_catalog([self.record(self.relative, self.source)])

    @staticmethod
    def record(relative, source):
        return {'path': relative, 'size_bytes': source.stat().st_size,
                'sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                'source_url': source.as_uri()}

    def write_catalog(self, records):
        self.catalog.write_text(json.dumps({'schema': 'mfm_bundled_assets_v1',
                                            'files': records, 'source_validation': {}}))

    def get(self, relative=None):
        return bundled_assets.get(self.config, self.root, relative or self.relative, required=True)

    def test_inventory_uses_packaged_catalog_without_creating_assets(self):
        base, records = bundled_assets.inventory(self.config, self.root)
        self.assertEqual(base, self.base)
        self.assertIn(self.relative, records)
        self.assertFalse(base.exists())

    def test_missing_asset_downloads_once_and_is_reused(self):
        with patch.object(asset_restore, 'fetch', wraps=asset_restore.fetch) as transfer:
            target = self.get()
            self.assertEqual(target.read_bytes(), self.source.read_bytes())
            self.assertEqual(self.get(), target)
        self.assertEqual(transfer.call_count, 1)

    def test_existing_valid_asset_does_not_transfer(self):
        target = self.base / self.relative
        target.parent.mkdir(parents=True)
        target.write_bytes(self.source.read_bytes())
        with patch.object(asset_restore, 'fetch', side_effect=AssertionError('중복 수신')):
            self.assertEqual(self.get(), target)

    def test_existing_corrupt_asset_is_preserved(self):
        target = self.base / self.relative
        target.parent.mkdir(parents=True)
        target.write_bytes(b'x' * self.source.stat().st_size)
        with patch.object(asset_restore, 'fetch') as transfer:
            with self.assertRaisesRegex(ValueError, 'SHA256 불일치'):
                self.get()
        transfer.assert_not_called()
        self.assertEqual(target.read_bytes(), b'x' * self.source.stat().st_size)

    def test_missing_fasta_and_index_are_prepared_in_dependency_order(self):
        fasta = self.root / 'source.fa'
        fasta.write_text('>chr1\nACGT\n')
        pysam.faidx(str(fasta))
        index = Path(str(fasta) + '.fai')
        relative = 'references/tiny.fa'
        self.write_catalog([self.record(relative, fasta), self.record(relative + '.fai', index)])
        target = self.get(relative + '.fai')
        self.assertEqual(target.read_bytes(), index.read_bytes())
        self.assertEqual((self.base / relative).read_bytes(), fasta.read_bytes())

    def test_shared_missing_asset_has_one_transfer_across_workers(self):
        with patch.object(asset_restore, 'fetch', wraps=asset_restore.fetch) as transfer:
            with ThreadPoolExecutor(max_workers=2) as pool:
                results = list(pool.map(lambda _: self.get(), range(2)))
        self.assertEqual(results[0], results[1])
        self.assertEqual(transfer.call_count, 1)

    def test_download_plan_does_not_treat_missing_asset_as_present(self):
        row = {'sample_id': 'S', 'read_type': 'short_read', 'platform': 'illumina',
               'alignment_path': '', 'rna_path': ''}
        resource = {'path': self.relative, 'url': self.source.as_uri(),
                    'kind': 'illumina_panel_vcf', 'shared': True, 'size_bytes': ''}
        config = dict(self.config, state_root=str(self.root / 'state'))
        with patch.object(download_only, 'sample_resources', return_value=[resource]):
            plan = download_only.plan(config, [row], self.root)
        self.assertEqual(plan['bundled_resources'], 0)
        self.assertEqual(plan['asset_restore_resources'], 1)
        self.assertEqual(plan['estimated_download_bytes'], self.source.stat().st_size)
        self.assertFalse(self.base.exists())

    def test_restore_dry_run_does_not_download(self):
        self.base.mkdir()
        (self.base / 'manifest.json').write_text(self.catalog.read_text())
        output = io.StringIO()
        with patch.object(asset_restore, 'fetch', side_effect=AssertionError('dry-run 수신')):
            with redirect_stdout(output):
                result = asset_restore.main(['--assets-root', str(self.base), '--group', 'panel', '--dry-run'])
        self.assertEqual(result, 0)
        self.assertEqual(json.loads(output.getvalue())['selected_files'], 1)
        self.assertFalse((self.base / self.relative).exists())

    def test_module_cli_runs_dry_run_without_downloading(self):
        self.base.mkdir()
        (self.base / 'manifest.json').write_text(self.catalog.read_text())
        module_root = Path(asset_restore.__file__).resolve().parents[1]
        source_root = Path(planner.__file__).resolve().parents[1]
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(
            [str(module_root), str(source_root), str(source_root / 'src')]))
        result = subprocess.run([sys.executable, '-m', 'handoff.asset_restore',
                                 '--assets-root', str(self.base), '--group', 'panel', '--dry-run'],
                                env=env, text=True, capture_output=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)['selected_files'], 1)
        self.assertFalse((self.base / self.relative).exists())


if __name__ == '__main__':
    unittest.main()
