from pathlib import Path
import sys,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from integrate_wgs_rna.pipeline.checkpoints import run_step,file_identity


class StepTest(unittest.TestCase):
    def test_reuse_invalidation_and_failure(self):
        with tempfile.TemporaryDirectory() as raw:
            r=Path(raw);source=r/'input';out=r/'output';marker=r/'step.json';source.write_text('a');calls=[]
            def action():calls.append(1);out.write_text(source.read_text());return {'status':'pass'}
            fn=lambda:run_step(marker,{'revision':1},action,lambda:file_identity(out),inputs_identity=lambda:file_identity(source))
            fn();fn();self.assertEqual(len(calls),1)
            out.write_text('wrong');fn();self.assertEqual(len(calls),2)
            source.write_text('b');fn();self.assertEqual(len(calls),3)
            def failed():raise RuntimeError('fixture failure')
            with self.assertRaisesRegex(RuntimeError,'fixture'):
                run_step(marker,{'revision':2},failed,lambda:file_identity(out))
            fn();self.assertEqual(len(calls),4)

    def test_changed_inputs_during_action_and_failed_result_are_not_reused(self):
        with tempfile.TemporaryDirectory() as raw:
            r=Path(raw);source=r/'input';source.write_text('a');marker=r/'step.json'
            def changing():source.write_text('b');return {'status':'pass'}
            with self.assertRaisesRegex(ValueError,'입력 지문 변경'):
                run_step(marker,{},changing,lambda:[],inputs_identity=lambda:file_identity(source))
            calls=[]
            def failure():calls.append(1);return {'status':'fail','reason':'실물 검사 실패'}
            for _ in range(2):
                self.assertEqual(run_step(marker,{},failure,lambda:[])['status'],'fail')
            self.assertEqual(len(calls),2)


if __name__ == '__main__':unittest.main()
