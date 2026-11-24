"""CLI의 local 연결·Slurm 분리·dry-run 무변경 계약을 검증함."""
from contextlib import redirect_stdout
from pathlib import Path
import copy
import io
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from handoff import cli


class QueueCliTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.config={'data_root':str(self.root/'input'),'output_root':str(self.root/'output'),
                     'state_root':str(self.root/'state'),'memory_limit_bytes':10**12,'model_target_bytes':10**12,
                     'local_jobs':2,'sample_cpus':8,'slurm':{'cpus':8,'concurrency':2}}
        self.rows=[{'sample_id':'S','estimated_model_bytes':'10','rna_path':''}]

    def invoke(self, command, executor='local'):
        args=['handoff',command,'--effective-config',str(self.root/'effective.json'),
              '--selection',str(self.root/'selection.tsv'),'--executor',executor]
        with patch.object(sys,'argv',args),patch.object(cli,'load_config',return_value=copy.deepcopy(self.config)), \
             patch.object(cli,'selection',return_value=self.rows):
            return cli.main()

    def test_local_run은공용큐로연결(self):
        with patch('handoff.local_queue.run',return_value=37) as run:
            self.assertEqual(self.invoke('run'),37)
        self.assertEqual(run.call_count,1)
        self.assertFalse((self.root/'state').exists())

    def test_slurm은기존제출경로유지(self):
        with patch('handoff.local_queue.run') as run,patch.object(cli.slurm,'submit',return_value={'submitted':False}) as submit,redirect_stdout(io.StringIO()):
            self.assertEqual(self.invoke('run','slurm'),0)
        run.assert_not_called();submit.assert_called_once()

    def test_doctor도실제공유풀계획을표시(self):
        output=io.StringIO()
        with redirect_stdout(output):self.assertEqual(self.invoke('doctor'),0)
        result=json.loads(output.getvalue())
        self.assertGreater(result['local_queue']['cpu_pool'],0)
        self.assertEqual(result['local_queue']['ram_reservation_bytes'],900_000_000_000)
        self.assertFalse((self.root/'state').exists())

    def test_dry_run은입력출력상태를생성하지않음(self):
        output=io.StringIO()
        with patch.object(cli,'sample_resources',return_value=[]),redirect_stdout(output):
            self.assertEqual(self.invoke('dry-run'),0)
        result=json.loads(output.getvalue())
        self.assertEqual(result['local_queue']['mode'],'ready_first_local')
        self.assertTrue(all(not Path(self.config[k]).exists() for k in ['data_root','output_root','state_root']))

if __name__=='__main__':unittest.main()
