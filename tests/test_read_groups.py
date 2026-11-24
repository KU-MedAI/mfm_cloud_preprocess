"""분담 보존·RAM예외·기존 진입점 전달 계약을 검사함."""
from pathlib import Path
from unittest.mock import patch
import contextlib
import io
import json
import tempfile
import unittest
from handoff import read_groups


class Tests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.config=read_groups.cli.load_config(read_groups.ROOT/'handoff_config.example.json')
        self.config['state_root']=str(Path(self.temp.name)/'state')
        self.config['data_root']=str(Path(self.temp.name)/'data')
        self.config['output_root']=str(Path(self.temp.name)/'output')

    def test_TGIL205명11배치_중복누락없음(self):
        items=read_groups.plans('tgil','long_read','all',self.config)
        self.assertEqual(len(items),11);self.assertEqual(sum(x['total'] for x in items),205)
        ids=[r['sample_id'] for x in items for r in x['rows']]+[s for x in items for s in x['blocked']]
        self.assertEqual(len(ids),len(set(ids)))
        for item in items:
            self.assertTrue(all(r['assigned_lab']=='TGIL' and r['batch_id']==item['batch'] and r['read_type']=='long_read' for r in item['rows']))
        self.assertFalse(Path(self.config['state_root']).exists())

    def test_MedAI001_89명과큰입력명시분리(self):
        self.config['pipeline_env']['MFM_LONGREAD_MEMORY_MODEL']='legacy'
        items=read_groups.plans('medai','long_read','001',self.config)
        self.assertEqual(items[0]['total'],89);self.assertIn('HG02554',items[0]['blocked'])
        self.assertNotIn('HG02554',[r['sample_id'] for r in items[0]['rows']])

    def test_short와long_state분리_output공유(self):
        long=read_groups.plans('tgil','long_read','001',self.config)[0]
        short=read_groups.plans('tgil','short_read','001',self.config)[0]
        self.assertNotEqual(long['config']['state_root'],short['config']['state_root'])
        self.assertEqual(long['config']['output_root'],short['config']['output_root'])
        self.assertEqual(long['config']['pipeline_env']['MFM_CPU_POOL_DIR'],short['config']['pipeline_env']['MFM_CPU_POOL_DIR'])

    def test_없는배치거부(self):
        with self.assertRaises(ValueError):read_groups.plans('tgil','long_read','999',self.config)

    def test_기존셸진입점으로run전달(self):
        self.config['pipeline_env']['MFM_LONGREAD_MEMORY_MODEL']='legacy'
        items=read_groups.plans('medai','long_read','001',self.config)
        with patch('handoff.read_groups.subprocess.run') as run,contextlib.redirect_stdout(io.StringIO()):
            run.return_value.returncode=0
            self.assertEqual(read_groups.run(items,self.config),1)
        command=run.call_args.args[0]
        self.assertEqual(command[:3],['bash',str(read_groups.ROOT/'run_pipeline.sh'),'run'])
        self.assertIn('--selection',command)
        summary=Path(self.config['state_root'])/'read_groups/medai_long_read_summary.json'
        self.assertEqual(json.loads(summary.read_text())['status'],'complete_with_pending')

    def test_실행오류후다음배치안함(self):
        items=read_groups.plans('tgil','short_read','all',self.config)
        with patch('handoff.read_groups.subprocess.run') as run:
            run.return_value.returncode=7
            self.assertEqual(read_groups.run(items,self.config),7)
            self.assertEqual(run.call_count,1)


if __name__=='__main__':unittest.main()
