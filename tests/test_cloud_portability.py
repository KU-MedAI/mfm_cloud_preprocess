"""클라우드 경로 격리·실제 RAM 한도·코드 배포 계약 회귀검사임."""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from handoff.storage_paths import storage_defaults
from handoff.resources import memory_budget
from tools.check_cloud_package import package_file, verify_package, export_package


class CloudPortabilityTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='mfm_cloud_test_')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_parent_storage_catalog_never_redirects_cloud_paths(self):
        parent = self.root / 'shared'
        package = parent / 'process_code/cloud'
        package.mkdir(parents=True)
        (package / '09_personal_fasta.sh').touch()
        catalog = parent / 'metadata/storage/paths.json'
        catalog.parent.mkdir(parents=True)
        catalog.write_text(json.dumps({'paths': [{'분류': '현재 MedAI 원본·캐시', '실제경로': str(parent / 'production')}]}))
        config = {'data_root': './data', 'output_root': './output', 'state_root': './run_state'}
        with patch.dict('os.environ', {'MFM_DATA_ROOT': str(parent)}):
            self.assertEqual(storage_defaults(config, package), config)
        self.assertFalse((parent / 'production').exists())

    def test_cloud_memory_budget_is_capped_by_actual_vm_allocation(self):
        with patch.dict('os.environ', {}, clear=True), patch('integrate_wgs_rna.memory_policy.allocation_memory_limit', return_value=64 << 30):
            self.assertEqual(memory_budget({'memory_limit_bytes': 1650000000000}), 64 << 30)
            self.assertEqual(memory_budget({'memory_limit_bytes': 32 << 30}), 32 << 30)

    def test_inventory_rejects_escape_and_symlinks(self):
        root = self.root / 'package'
        root.mkdir()
        outside = self.root / 'outside'
        outside.write_text('private')
        (root / 'link').symlink_to(outside)
        for value in ['../outside', str(outside), 'link']:
            with self.subTest(value=value), self.assertRaises(ValueError):
                package_file(root, value)

    def test_export_excludes_unlisted_raw_and_refuses_existing_directory(self):
        root = self.root / 'source'
        root.mkdir()
        payload = root / 'code.py'
        payload.write_text('print("ok")\n')
        raw = root / 'data/raw.bam'
        raw.parent.mkdir()
        raw.write_bytes(b'raw must stay here')
        records = [{'path': payload.name, 'size_bytes': payload.stat().st_size,
                    'sha256': hashlib.sha256(payload.read_bytes()).hexdigest()}]
        (root / 'CLOUD_MANIFEST.json').write_text(json.dumps({'schema': 'mfm_cloud_package_v1', 'files': records}))
        destination = self.root / 'exported'
        with patch('tools.check_cloud_package.REQUIRED', {payload.name}):
            self.assertEqual(export_package(root, destination)['status'], 'pass')
            self.assertFalse((destination / 'data').exists())
            with self.assertRaisesRegex(ValueError, '이미 존재'):
                export_package(root, destination)
            (destination / 'code.py').write_text('changed')
            self.assertEqual(verify_package(destination)[0]['status'], 'fail')

    def test_preflight_needs_only_python_standard_library(self):
        package = Path(__file__).resolve().parents[1]
        env = {'PATH': os.defpath, 'HOME': str(self.root), 'PYTHONDONTWRITEBYTECODE': '1',
               'PYTHONPATH': os.pathsep.join(str(package / p) for p in ['.', 'src', 'vendor/MFM-DNARNAFM'])}
        for command in ['doctor', 'dry-run', 'download']:
            with self.subTest(command=command):
                args = [sys.executable, '-S', '-m', 'handoff.cli', command,
                        '--config', str(package / 'cloud_config.pileup.example.json')]
                if command != 'doctor':
                    args += ['--lab', 'medai', '--batch', '001', '--sample', 'NA12878']
                if command == 'download':
                    args += ['--dry-run']
                result = subprocess.run(args, cwd=package, env=env, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIsInstance(json.loads(result.stdout), dict)
        self.assertFalse((self.root / '.conda_miniforge').exists())

    def test_deployed_profiles_have_no_research_server_path_dependency(self):
        package = Path(__file__).resolve().parents[1]
        profiles = [package / 'handoff_config.example.json', package / 'cloud_config.pileup.example.json']
        profiles += sorted((package / 'manifests/클라우드_분담/설정').glob('*.json'))
        self.assertGreaterEqual(len(profiles), 6)
        for path in profiles:
            with self.subTest(profile=path.name):
                config = json.loads(path.read_text())
                for key in ['data_root', 'output_root', 'state_root', 'bundled_assets_root']:
                    self.assertIn(package, (path.parent / config[key]).resolve().parents)
                policy = config.get('pipeline_env', {})
                if policy.get('BLACKLIST_BED'):
                    self.assertTrue((package / policy['BLACKLIST_BED']).is_file())
                if policy.get('MFM_CPU_POOL_DIR'):
                    self.assertIn(package, (package / policy['MFM_CPU_POOL_DIR']).resolve().parents)
                self.assertEqual(config['slurm']['cpus'], config['sample_cpus'])

    def test_invalid_inventory_records_report_failure(self):
        root = self.root / 'invalid'
        root.mkdir()
        for records in ['wrong', [None], [{'path': 'x', 'size_bytes': 'wrong', 'sha256': 'x'}]]:
            (root / 'CLOUD_MANIFEST.json').write_text(json.dumps({'schema': 'mfm_cloud_package_v1', 'files': records}))
            self.assertEqual(verify_package(root)[0]['status'], 'fail')
