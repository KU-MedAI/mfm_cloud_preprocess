"""합산 예약·동시 실행·캐시 재사용과 실패 시 완료 장부 생성을 확인한다."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from threading import Barrier
import json
import os
import tempfile
import unittest
from integrate_wgs_rna.gpu_prepare import prepare_short
from integrate_wgs_rna.memory_policy import GIB


class PreparationTests(unittest.TestCase):
    def test_GPU대기중_실행중인생산자_CPU를반환하지않음(self):
        from integrate_wgs_rna.gpu_stage_resources import StageBudget
        env={'MFM_GPU_STAGE_CPUS':'68','MFM_GPU_COMPARE_UNRESERVED':'1',
             'RUN_ID':'gpu_compare_live_input','MFM_WGS_BACKEND':'parabricks'}
        with patch.dict(os.environ,env,clear=True):
            budget=StageBudget();closed=[];lease=SimpleNamespace(close=lambda:closed.append(True))
            with patch('integrate_wgs_rna.gpu_pool.acquire',return_value=lease),\
                 patch.object(budget,'available_cpus',side_effect=[32,68]),patch('time.sleep'):
                with budget.gpu(1024,purpose='test',producer_cpus=32):pass
            self.assertEqual(budget.events[0]['phase'],'gpu_wait_input_active')
            self.assertEqual(budget.events[0]['cpus'],32)
            self.assertEqual(budget.events[-1]['cpus'],68)
            self.assertEqual(len(closed),2)

    def test_같은예약안에서실제동시진입(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ref=root/'ref.fa';ref.write_text('>chr1\nACGT\n')
            events=[];barrier=Barrier(2,timeout=5)
            budget=SimpleNamespace(maximum_cpus=16,available_cpus=lambda:16,memory=128*GIB,
                phase=lambda *a,**k:events.append((a,k)),
                try_phase=lambda *a,**k:events.append((a,k)) or True)
            def index(reference,work):
                barrier.wait();return reference
            def inputs(source,key,work,paired,share,producer,expected):
                barrier.wait();self.assertEqual(share.maximum_cpus,15)
                share.wait()
                share.phase('input',15,4*GIB)
                with self.assertRaises(ValueError):share.phase('over',16,4*GIB)
                with self.assertRaises(ValueError):share.phase('over',15,5*GIB)
                return [root/'r1.gz',root/'r2.gz'],2,False
            with patch('integrate_wgs_rna.gpu_prepare.cached_inputs',side_effect=inputs):
                result=prepare_short(ref,root,'source','key',root,budget,3,100,index)
            self.assertEqual(len(events),4)
            self.assertEqual(events[0][0][0],'bwa_index_prepare')
            self.assertEqual(events[1][0][0],'input_cache_wait_index_active')
            self.assertEqual(events[1][0][1],1)
            self.assertEqual(events[2][0][0],'bwa_index_and_input')
            self.assertEqual(events[3][0][0],'input_cache_wait_index_active')
            self.assertEqual(events[3][0][1],1)
            plan=result[-1]
            self.assertEqual(plan['index_cpus']+plan['input_cpus'],16)
            self.assertEqual(events[0][0][2],plan['reserved_memory_bytes'])
            self.assertEqual(events[1][0][2],events[0][0][2])
            self.assertGreater(plan['expected_disk_peak_bytes'],8*GIB)

    def test_부족한CPU_RAM은기존순차경로로반환(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ref=root/'r';ref.write_bytes(b'A')
            for cpus,memory in [(4,128*GIB),(16,GIB)]:
                budget=SimpleNamespace(maximum_cpus=cpus,available_cpus=lambda:cpus,memory=memory)
                self.assertIsNone(prepare_short(ref,root,'s','k',root,budget,3,100,None))

    def test_main은동시준비결과를정렬과메타데이터에사용(self):
        from contextlib import nullcontext
        from integrate_wgs_rna import gpu_wgs
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ,{},clear=True):
            root=Path(tmp);ref=root/'ref';ref.write_text('>chr1\nACGT\n')
            fastqs=[root/'r1.gz',root/'r2.gz']
            for f in fastqs:f.write_bytes(b'data')
            plan={'index_seconds':2,'input_seconds':3,'elapsed_seconds':3}
            budget=SimpleNamespace(cpus=16,maximum_cpus=16,events=[],memory=128*GIB,
                phase=lambda *a,**k:None,gpu=lambda *a,**k:nullcontext(None),output=lambda:None)
            def run(command,*args,**kwargs):
                Path(command[command.index('--out-bam')+1]).write_bytes(b'BAM')
                self.assertEqual(command[command.index('--in-fq')+1:command.index('--in-fq')+3],list(map(str,fastqs)))
                return {}
            runtime=SimpleNamespace(validate=lambda:None,run=run)
            args=['--reference',str(ref),'--platform','illumina','--preset','sr','--sample','S',
                  '--read-group','@RG\\tID:S','--work',str(root),'--metadata',str(root/'meta.json'),
                  '--threads','16','--input-source',str(root/'source.sh'),
                  '--input-cache-key',json.dumps({'alignment':['unused',100]})]
            with patch.object(gpu_wgs,'Runtime',return_value=runtime),patch.object(gpu_wgs,'StageBudget',return_value=budget),\
                 patch('integrate_wgs_rna.gpu_bwa_index.lookup',return_value=None),\
                 patch('integrate_wgs_rna.gpu_prepare.prepare_short',return_value=(ref,fastqs,2,False,plan)),\
                 patch.object(gpu_wgs,'restore_bam',return_value={}),\
                 patch.object(gpu_wgs,'bwa_reference',side_effect=AssertionError('두번 준비하면 안 됨')):
                self.assertEqual(gpu_wgs.main(args),0)
            meta=json.loads((root/'meta.json').read_text())
            self.assertEqual(meta['input_preparation_overlap'],plan)
            self.assertFalse(meta['input_streamed'])
            self.assertEqual(meta['input_records'],2)


if __name__=='__main__':unittest.main()
