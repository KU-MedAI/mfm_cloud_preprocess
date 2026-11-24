"""설치 전 계획과 배포 Slurm 설정이 실제 자원 계약을 만족하는지 확인함."""
from pathlib import Path
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


class BootstrapContractTests(unittest.TestCase):
    def test_preflight_works_without_installed_python_packages(self):
        commands = [['doctor']]
        commands += [['dry-run', '--target', name] for name in
                     ['azure_poc', 'azure_main', 'aws_main', 'server']]
        commands += [['download', '--target', 'azure_poc', '--sample', 'HG00308', '--dry-run']]
        with tempfile.TemporaryDirectory() as home:
            env = {'PATH': os.defpath, 'HOME': home, 'PYTHONDONTWRITEBYTECODE': '1',
                   'PYTHONPATH': os.pathsep.join(str(ROOT / p) for p in
                                               ['.', 'src', 'vendor/MFM-DNARNAFM'])}
            for command in commands:
                with self.subTest(command=command):
                    result = subprocess.run([sys.executable, '-S', '-m', 'handoff.cli', *command],
                                            cwd=ROOT, env=env, capture_output=True, text=True, timeout=30)
                    self.assertEqual(result.returncode, 0, result.stderr[-3000:])
                    self.assertIsInstance(json.loads(result.stdout), dict)
            self.assertFalse((Path(home) / '.conda_miniforge').exists())

    def test_slurm_profiles_fit_two_star_haps_and_parent_with_margin(self):
        from handoff.slurm import commands
        from integrate_wgs_rna.memory_policy import GIB, usable_memory
        from integrate_wgs_rna.rna_memory import estimate_star_memory, apply_star_memory
        profiles = [ROOT / 'handoff_config.example.json']
        profiles += sorted((ROOT / 'manifests/클라우드_분담/설정').glob('*.json'))
        extra = ROOT / 'cloud_config.pileup.example.json'
        if extra.exists():
            profiles.append(extra)
        for profile in profiles:
            with self.subTest(profile=profile), tempfile.TemporaryDirectory() as raw:
                config = json.loads(profile.read_text())
                request = config['slurm']['memory']
                self.assertTrue(request.endswith('G'))
                gib = int(request[:-1])
                plans = [estimate_star_memory(mode, config['slurm']['cpus'],
                                              reference_bases=4_000_000_000)
                         for mode in ['genomeGenerate', 'alignReads']]
                plan = max(plans, key=lambda p: p['memory_bytes'])
                pool = Path(raw); (pool / 'leases').mkdir()
                (pool / 'policy.json').write_text(json.dumps({'memory_limit_bytes': gib * GIB}))
                (pool / 'leases/parent.json').write_text(json.dumps({'memory_bytes': 4 * GIB}))
                env = {'SLURM_JOB_ID': 'test-only', 'SLURM_MEM_PER_NODE': str(gib * 1024),
                       'SLURM_CPUS_PER_TASK': str(config['slurm']['cpus'])}
                with patch.dict(os.environ, env, clear=True):
                    required = 2 * plan['memory_bytes'] + 4 * GIB
                    self.assertGreaterEqual(usable_memory(pool), required,
                                            '두 HAP STAR와 부모 예약에 10% 여유가 필요함')
                    for hap in [1, 2]:
                        lease = pool / f'leases/hap{hap}.json'
                        lease.write_text(json.dumps({'branch': 'RNA', 'memory_bytes': 4 * GIB}))
                        self.assertEqual(apply_star_memory(lease, plan)['memory_bytes'], plan['memory_bytes'])
                config.update(state_root=str(pool), slurm=dict(config['slurm']))
                array = commands(config, ROOT, pool / 'config.json', pool / 'selection.tsv', 1)[1]
                self.assertIn('--mem=' + request, array)


if __name__ == '__main__':
    unittest.main()
