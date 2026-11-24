"""완료 정렬 복구의 성공 증거와 변경/실패 거부를 검증함."""
import json
import os
from pathlib import Path
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from integrate_wgs_rna.realign.recover_output import completed_output_proof


class RecoveryGuardTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.input = self.root/'input.bam'
        self.target = self.root/'target.fa'
        self.output = self.root/'out.cram'
        self.script = self.root/'TEST.HAP1.realign.sh'
        self.resource = self.root/'TEST.HAP1.resource.json'
        for path in [self.input, self.target, self.output]:
            path.write_text('guard fixture only')
        self.script.write_text('\n'.join([
            f'INPUT={self.input}', f'TARGET={self.target}',
            f'OUTPUT={self.output}', f'PART={self.output}.part',
            'minimap2 -ax map-hifi', 'mv "$PART" "$OUTPUT"']))
        self.started = time.time()-100
        for path in [self.input, self.target, self.script]:
            os.utime(path, (self.started-10, self.started-10))
        os.utime(self.output, (self.started+10, self.started+10))
        self.usage = dict(exit_code=0, command=['bash', str(self.script)],
                          started_at_unix=self.started, elapsed_sec=11,
                          user_cpu_sec=40, system_cpu_sec=2)
        self.resource.write_text(json.dumps(self.usage))
        self.config = SimpleNamespace(overwrite=False, dry_run=False, platform='hifi')
        self.env = patch.dict(os.environ, MFM_REALIGN_RESUME_COMPLETED_OUTPUT='1', MFM_WGS_BACKEND='cpu')
        self.env.start()
        self.addCleanup(self.env.stop)

    def proof(self):
        return completed_output_proof(self.config, 'TEST', 1, self.input,
                                      self.target, self.output, self.root)

    def test_successful_completed_output(self):
        proof = self.proof()
        self.assertEqual(proof['status'], 'verified_completed_alignment')
        self.assertEqual(proof['original_cpu_seconds'], 42)

    def test_default_does_not_reuse(self):
        with patch.dict(os.environ, MFM_REALIGN_RESUME_COMPLETED_OUTPUT='0'):
            self.assertIsNone(self.proof())

    def test_missing_output_uses_normal_alignment(self):
        self.output.unlink()
        self.assertIsNone(self.proof())

    def test_failed_command_rejected(self):
        self.usage['exit_code'] = 1
        self.resource.write_text(json.dumps(self.usage))
        with self.assertRaises(ValueError):
            self.proof()

    def test_changed_input_reference_or_script_rejected(self):
        for path in [self.input, self.target, self.script]:
            with self.subTest(path=path.name):
                os.utime(path, (self.started+20, self.started+20))
                with self.assertRaises(ValueError):
                    self.proof()
                os.utime(path, (self.started-10, self.started-10))

    def test_wrong_input_path_rejected(self):
        self.script.write_text(self.script.read_text().replace(f'INPUT={self.input}', 'INPUT=/other.bam'))
        os.utime(self.script, (self.started-10, self.started-10))
        with self.assertRaises(ValueError):
            self.proof()

    def test_partial_output_rejected(self):
        Path(str(self.output)+'.part').touch()
        with self.assertRaises(ValueError):
            self.proof()

    def test_overwrite_and_dry_run_rejected(self):
        for key in ['overwrite', 'dry_run']:
            with self.subTest(key=key):
                setattr(self.config, key, True)
                with self.assertRaises(ValueError):
                    self.proof()
                setattr(self.config, key, False)


if __name__ == '__main__':
    unittest.main()
