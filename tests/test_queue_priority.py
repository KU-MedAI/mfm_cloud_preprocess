"""실제 큐 Harness로 drain·long우선·버퍼·재개 식별 보존을 검사함."""
from pathlib import Path
import tempfile
import unittest
from test_ready_queue import Harness, Process


class Tests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)

    def harness(self, *names):
        h = Harness(self.temp.name, names)
        for i, row in enumerate(h.rows):
            row.update(read_type='long_read' if row['sample_id'].startswith('long') else 'short_read',
                       known_input_total_bytes=str((i+1)*100))
        h.queue.PRIORITY_POLICY = {'mode': 'long_read_first'}
        return h

    def test_기존6명모두끝난뒤계산_다운로드는선행(self):
        h=self.harness('short_old','short_next','long_first');h.success('short_old');h.success('short_next')
        h.queue.PRIORITY_POLICY['drain_samples']=['short_old']
        s=h.scheduler();s.active['short_old']=(s.entries[0],Process())
        s.step();self.assertEqual(set(s.active),{'short_old'});self.assertIn('long_first',s.downloads)
        h.success('long_first');s.downloads['long_first'][1].result=0
        s.step();self.assertEqual(set(s.active),{'short_old'})
        h.done.add('short_old');s.active['short_old'][1].result=0;s.step()
        self.assertEqual(set(s.active),{'long_first'});self.assertNotIn('short_next',s.active)

    def test_short준비버퍼가가득차도long다운로드진행(self):
        h=self.harness(*[f'short{i}' for i in range(12)],'long_next')
        for name in [f'short{i}' for i in range(12)]:h.success(name)
        s=h.scheduler();s.step();self.assertIn('long_next',s.downloads);self.assertFalse(s.active)

    def test_long완료후short재개_큰long보류도기록(self):
        h=self.harness('short','long_small','long_big')
        for row in h.rows:h.success(row['sample_id'])
        h.queue.PRIORITY_POLICY['blocked_samples']={'long_big':'RAM 초과'}
        s=h.scheduler();s.step();self.assertEqual(set(s.active),{'long_small'})
        self.assertNotIn('long_big',s.downloads)
        h.done.add('long_small');s.active['long_small'][1].result=0;s.step();self.assertEqual(set(s.active),{'short'})
        event=[d for n,d in h.events if n=='sample_ready_queue'][-1]
        self.assertIn('long_big',event['deferred_samples']);self.assertEqual(event['priority_policy']['blocked_samples'],{'long_big':'RAM 초과'})

    def test_다운로드중short와작업인덱스유지(self):
        h=self.harness('short_first','long_next');s=h.scheduler()
        before=list(s.entries);s.downloads['short_first']=(s.entries[0],Process());s.step()
        self.assertIn('short_first',s.downloads);self.assertIn('long_next',s.downloads);self.assertEqual(s.entries,before)
        started=[(k,n) for k,n,_ in h.started];self.assertNotIn(('worker','short_first'),started)

    def test_long실패격리후short진행(self):
        h=self.harness('short','long_bad');h.success('short');s=h.scheduler()
        s.ledger['downloads']['long_bad']={'attempts':3,'deferred':True,'next_retry_at':0}
        s.step();self.assertIn('short',s.active)

    def test_고정모드_실행중정책으로local_jobs상한을늘리거나해제못함(self):
        for requested in [None,4,64]:
            with self.subTest(requested=requested):
                h=self.harness('short1','short2','short3')
                h.queue.config={'cpu_allocation_mode':'fixed_sample','local_jobs':2}
                h.queue.PRIORITY_POLICY={'mode':'ready_first','max_active_samples':requested}
                s=h.scheduler()
                self.assertEqual(2,h.queue.MAX_ACTIVE_SAMPLES)
                for r in h.rows:h.success(r['sample_id'])
                s.step();self.assertEqual(2,len(s.active))

    def test_고정모드_더작은정책상한과실행중변경보존(self):
        h=self.harness('short1','short2','short3')
        h.queue.config={'cpu_allocation_mode':'fixed_sample','local_jobs':2}
        h.queue.PRIORITY_POLICY={'mode':'ready_first','max_active_samples':1}
        s=h.scheduler();self.assertEqual(1,h.queue.MAX_ACTIVE_SAMPLES)
        h.save(h.base/'우선처리정책.json',{'mode':'ready_first','max_active_samples':None})
        s.priority.refresh();self.assertEqual(2,h.queue.MAX_ACTIVE_SAMPLES)

    def test_동적모드의명시적동시수해제는유지(self):
        h=self.harness('short')
        h.queue.config={'cpu_allocation_mode':'dynamic','local_jobs':2}
        h.queue.PRIORITY_POLICY={'mode':'ready_first','max_active_samples':None}
        h.scheduler();self.assertIsNone(h.queue.MAX_ACTIVE_SAMPLES)

    def test_설정오류거부(self):
        h=self.harness('short');h.queue.PRIORITY_POLICY={'mode':'invalid'}
        with self.assertRaises(ValueError):h.scheduler()

    def test_설정파일원자교체후정책갱신(self):
        h=self.harness('short','long');s=h.scheduler()
        h.save(h.base/'우선처리정책.json',{'mode':'ready_first'})
        h.success('short');s.step();self.assertIn('short',s.active)

    def test_long우선입력작은순서(self):
        h=self.harness('long_big','long_small');h.rows[0]['known_input_total_bytes']='2000'
        s=h.scheduler();s.step();names=[n for k,n,_ in h.started if k=='download-worker']
        self.assertEqual(names,['long_small','long_big'])

    def test_실패한기존drain작업은전체큐막지않음(self):
        h=self.harness('short_old','long');h.success('long')
        h.queue.PRIORITY_POLICY['drain_samples']=['short_old'];s=h.scheduler()
        s.ledger['workers']['short_old']={'reason':'실패'};s.step();self.assertIn('long',s.active)

    def test_long전용은_마지막long끝나도short진입금지(self):
        h=self.harness('short','long');h.queue.PRIORITY_POLICY={'mode':'long_read_only'}
        h.success('short');h.success('long');s=h.scheduler();s.step()
        self.assertEqual(set(s.active),{'long'})
        h.done.add('long');s.active['long'][1].result=0
        self.assertEqual(s.step(),1);self.assertFalse(s.active);self.assertFalse(s.downloads)
        self.assertNotIn('short',s.ledger['workers'])
        self.assertTrue(s.priority.paused('short'))

    def test_long전용_자원보류만남아도short다운로드안함(self):
        h=self.harness('short','long');h.queue.PRIORITY_POLICY={'mode':'long_read_only','blocked_samples':{'long':'RAM초과'}}
        s=h.scheduler();self.assertEqual(s.step(),1);self.assertFalse(h.started)
        status=[x for n,x in h.events if n=='sample_ready_queue'][-1]
        self.assertEqual(status['deferred_samples'],['long'])
        self.assertEqual(status['priority_policy']['paused_short_samples'],['short'])

    def test_long전용에서도기존short작업6명은보존(self):
        h=self.harness('short_old','short_next','long');h.queue.PRIORITY_POLICY={'mode':'long_read_only','drain_samples':['short_old']}
        h.success('short_old');h.success('short_next');s=h.scheduler()
        process=Process();s.active['short_old']=(s.entries[0],process);s.step()
        self.assertEqual(set(s.active),{'short_old'});self.assertIn('long',s.downloads)
        self.assertFalse(s.priority.paused('short_old'));self.assertTrue(s.priority.paused('short_next'))


    def test_기존작업보존하면서빈자원에long즉시투입(self):
        h=self.harness('short_old','short_next','long');h.queue.PRIORITY_POLICY={'mode':'long_read_only','drain_samples':['short_old'],'wait_for_drain':False}
        h.success('short_next');h.success('long');s=h.scheduler();s.active['short_old']=(s.entries[0],Process())
        s.step();self.assertEqual(set(s.active),{'short_old','long'});self.assertNotIn('short_next',s.active)
        self.assertEqual(s.priority.status()['drain_waiting'],[])
        self.assertEqual(s.priority.status()['preserved_active_samples'],['short_old'])





    def test_long_pair는큰입력이어도WGS복구보다우선(self):
        h=self.harness('long_wgs','long_pair');h.rows[0]['data_group']='wgs';h.rows[1]['data_group']='wgs_rna_pair';h.rows[1]['known_input_total_bytes']='999999'
        s=h.scheduler();s.recoveries['long_wgs']={'priority':-100}
        for worker in [False,True]:
            self.assertEqual(s.priority.ordered(s.entries,worker=worker)[0][2]['sample_id'],'long_pair')

    def test_pair_first설정형식검사(self):
        h=self.harness('long');h.queue.PRIORITY_POLICY['pair_first']='yes'
        with self.assertRaises(ValueError):h.scheduler()


    def test_동시상한실시간축소_진행6명보존_3명후1명투입(self):
        h=self.harness(*[f'long{i}' for i in range(8)])
        h.queue.MAX_ACTIVE_SAMPLES=None
        for row in h.rows:h.success(row['sample_id'])
        s=h.scheduler()
        for entry in s.entries[:6]:s.active[entry[2]['sample_id']]=(entry,Process())
        preserved={name:process.pid for name,(_,process) in s.active.items()}
        h.save(h.base/'우선처리정책.json',{'mode':'long_read_only','max_active_samples':4})
        s.step()
        self.assertEqual(len(s.active),6)
        self.assertFalse([x for x in h.started if x[0]=='worker'])
        self.assertEqual(s.priority.status()['max_active_samples'],4)
        for name in ['long0','long1','long2']:
            h.done.add(name);s.active[name][1].result=0
        s.step()
        self.assertEqual(set(s.active),{'long3','long4','long5','long6'})
        for name in ['long3','long4','long5']:
            self.assertEqual(s.active[name][1].pid,preserved[name])

    def test_동시상한실시간변경_해제시시작설정복원(self):
        h=self.harness('long0','long1','long2');h.queue.MAX_ACTIVE_SAMPLES=None
        for row in h.rows:h.success(row['sample_id'])
        h.save(h.base/'우선처리정책.json',{'mode':'long_read_only','max_active_samples':1})
        s=h.scheduler();s.step();self.assertEqual(len(s.active),1)
        h.save(h.base/'우선처리정책.json',{'mode':'long_read_only','max_active_samples':2})
        s.step();self.assertEqual(len(s.active),2)
        h.save(h.base/'우선처리정책.json',{'mode':'long_read_only'})
        s.step();self.assertEqual(len(s.active),3)
        self.assertIsNone(h.queue.MAX_ACTIVE_SAMPLES)

    def test_동시상한은다운로드를막지않음(self):
        h=self.harness('long0','long1','long2');h.success('long0');h.success('long1')
        h.queue.DOWNLOAD_TARGET=4
        h.queue.PRIORITY_POLICY['max_active_samples']=1
        s=h.scheduler();s.step()
        self.assertEqual(len(s.active),1)
        self.assertIn('long2',s.downloads)

    def test_동시상한형식검사(self):
        for invalid in [0,-1,True,'4',4.5]:
            with self.subTest(value=invalid):
                h=self.harness('long');h.queue.PRIORITY_POLICY['max_active_samples']=invalid
                with self.assertRaises(ValueError):h.scheduler()
