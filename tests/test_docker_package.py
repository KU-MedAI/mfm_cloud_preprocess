"""Docker 진입점의 인수·캐시·중단 처리와 실제 파이프라인 연결을 확인함."""
from pathlib import Path
import json
import os
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class DockerPackageTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='mfm_docker_test_')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        package = self.root / 'package'
        (package / 'docker').mkdir(parents=True)
        self.entry = package / 'docker/entrypoint.sh'
        shutil.copy2(ROOT / 'docker/entrypoint.sh', self.entry)
        pipeline = package / 'run_pipeline.sh'
        pipeline.write_text('#!/usr/bin/env bash\nexec "$TEST_PYTHON" -c \'import json,os,signal,sys; '
                            'print(json.dumps({"args":sys.argv[1:],"home":os.environ["HOME"]}),flush=True); '
                            'signal.pause() if "wait" in sys.argv else None\' "$@"\n')
        self.seed = self.root / 'seed'; self.seed.mkdir()
        self.cache = self.root / 'cache/native'
        self.env = dict(os.environ, HOME=str(self.root / 'home'), TMPDIR=str(self.root / 'tmp'),
                        XDG_CACHE_HOME=str(self.root / 'cache'), MFM_NATIVE_CACHE_ROOT=str(self.cache),
                        MFM_IMAGE_NATIVE_CACHE=str(self.seed), TEST_PYTHON=sys.executable,
                        PYTHONDONTWRITEBYTECODE='1')

    def invoke(self, *args, env=None):
        return subprocess.run(['bash', str(self.entry), *args], env=env or self.env,
                              capture_output=True, text=True, timeout=20)

    def test_entrypoint_forwards_arguments_and_defaults_to_help(self):
        for args, expected in [((), ['help']), (('run', '--sample', 'space sample'),
                                              ['run', '--sample', 'space sample'])]:
            with self.subTest(args=args):
                result = self.invoke(*args)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(json.loads(result.stdout)['args'], expected)
                self.assertTrue(self.cache.is_dir())
                self.assertTrue(Path(self.env['HOME']).is_dir())

    def test_native_seed_reuses_build_and_preserves_existing_cache(self):
        seed = self.seed / 'native_key'; seed.write_bytes(b'image native')
        (self.seed / 'native_key.lock').write_text('build lock')
        result = self.invoke('help')
        self.assertEqual(result.returncode, 0, result.stderr)
        target = self.cache / seed.name
        self.assertEqual(target.read_bytes(), b'image native')
        self.assertEqual(stat.S_IMODE(target.stat().st_mode), 0o700)
        self.assertEqual((self.cache / 'native_key.lock').stat().st_size, 0)
        target.write_bytes(b'existing native')
        self.assertEqual(self.invoke('help').returncode, 0)
        self.assertEqual(target.read_bytes(), b'existing native')

    def test_failed_seed_copy_does_not_publish_a_partial_binary(self):
        (self.seed / 'native_key').write_bytes(b'image native')
        # /tmp가 noexec인 서버에서도 실패한 복사를 실제로 모의함.
        env = dict(self.env)
        env['BASH_FUNC_cp%%'] = '() { printf partial > "${@: -1}"; return 1; }'
        result = self.invoke('help', env=env)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.cache / 'native_key').exists())
        self.assertFalse(list(self.cache.glob('.seed.*')))

    def test_exec_delivers_sigterm_to_the_pipeline(self):
        process = subprocess.Popen(['bash', str(self.entry), 'wait'], env=self.env,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            self.assertEqual(json.loads(process.stdout.readline())['args'], ['wait'])
            process.terminate()
            process.communicate(timeout=10)
            self.assertEqual(process.returncode, -signal.SIGTERM)
        finally:
            if process.poll() is None:
                process.kill(); process.communicate()

    def test_docker_config_preserves_biology_and_limits_fixed_cpu(self):
        original = json.loads((ROOT / 'cloud_config.pileup.example.json').read_text())
        docker = json.loads((ROOT / 'docker/config.pileup.example.json').read_text())
        policy = dict(docker['pipeline_env']); policy.pop('MFM_SCRATCH_ROOT')
        self.assertEqual(policy, original['pipeline_env'])
        for key in ['remove_owned_raw_inputs', 'cleanup_shared_fasta', 'acceleration']:
            self.assertEqual(docker[key], original[key])
        self.assertLessEqual(docker['sample_cpus'] * docker['local_jobs'], docker['cpu_pool_limit'])
        for key in ['data_root', 'output_root', 'state_root', 'bundled_assets_root']:
            self.assertTrue(docker[key].startswith('/data/'))

    def test_entrypoint_connects_to_real_preflight_without_installing(self):
        env = dict(self.env, PYTHONPATH=os.pathsep.join(str(ROOT / p) for p in
                                                      ['.', 'src', 'vendor/MFM-DNARNAFM']))
        result = subprocess.run(['bash', str(ROOT / 'docker/entrypoint.sh'), 'dry-run',
                                 '--target', 'azure_poc', '--sample', 'HG00308',
                                 '--config', str(ROOT / 'docker/config.pileup.example.json'),
                                 '--sample-cpus', '1', '--cpu-pool-limit', '1'],
                                cwd=ROOT, env=env, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads(result.stdout)
        self.assertEqual(report['sample_ids'], ['HG00308'])
        self.assertEqual(report['config']['data_root'], '/data/input')
        self.assertEqual(report['config']['sample_cpus'], 1)
        self.assertEqual(report['config']['pipeline_env']['MFM_LONG_READ_OUTPUT_MODE'], 'pileup_feature')
        self.assertFalse((self.root / '.conda_miniforge').exists())


if __name__ == '__main__':
    unittest.main()
