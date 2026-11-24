"""실제 다운로드 없이 샘플 실패 격리·재시도·후속 진행을 검증함."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/"src")]
from types import SimpleNamespace
import hashlib
import json
import tempfile
import unittest

from handoff.ready_queue import Scheduler, ExistingDownload


class Process:
    counter = 100
    def __init__(self, result=None):
        Process.counter += 1
        self.pid = Process.counter
        self.result = result
    def poll(self):
        return self.result
    def wait(self):
        return self.result


class Harness:
    def __init__(self, base, names):
        self.base = Path(base)
        self.clock = 1000
        self.rows = [{'sample_id': name} for name in names]
        self.work = [(i + 1, [row]) for i, row in enumerate(self.rows)]
        self.done = set()
        self.records = {}
        self.events = []
        self.started = []
        self.allowed = True
        self.locks = set()
        common = SimpleNamespace(
            save=self.save, files=self.files, completed_ids=lambda rows: set(self.done),
            verified=lambda row, state: row['sample_id'] in self.done,
            event=lambda status, **details: self.events.append((status, details)))
        self.queue = SimpleNamespace(
            common=common, BASE=self.base, MAX_ACTIVE_SAMPLES=4, state_root=lambda path: Path(json.loads(path.read_text())['state_root']),
            download_record=lambda state, wave, sample: self.records.get(sample, {}),
            download_state=lambda state, wave: state / 'download_only',
            sample_locked=lambda state, sample: (str(state), sample) in self.locks,
            resource_check=lambda *a, **kw: (self.allowed, '' if self.allowed else '디스크 부족'),
            prepare=lambda *a: None, start=self.start, worker_config=lambda context, sample, grant: context['config'],
            cpu_quota=lambda: 68, reserved_ram=lambda rows: len(rows) * 10)

    def save(self, path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))

    def files(self, number, wave):
        path = self.base / f'config{number}.json'
        self.save(path, {'state_root': str(self.base / f'state{number}'), 'sample_cpus': 68})
        return str(number), self.base / f'select{number}.tsv', path

    def start(self, kind, tag, selection, config, **kwargs):
        row = self.work[int(tag) - 1][1][0]
        process = Process(0 if kind == 'finish' else None)
        self.started.append((kind, row['sample_id'], process))
        return process

    def scheduler(self):
        s = Scheduler(self.queue, self.work, self.rows, clock=lambda: self.clock)
        s.initialize_download = lambda ctx: None
        return s

    def success(self, sample):
        self.records[sample] = {'status': 'success', 'files': []}


class Tests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)

    def harness(self, *names):
        return Harness(self.directory.name, names)

    def test_모델입력목표도달시_신규다운로드와계산중지(self):
        h=self.harness('done','next');h.done.add('done');h.success('next')
        h.queue.target_reached=lambda:True
        self.assertEqual(h.scheduler().step(),0)
        self.assertFalse(h.started)
        self.assertEqual(h.events[-1][0],'target_reached')

    def test_공용큐의명시다운로드와버퍼설정존중(self):
        h=self.harness('a','b','c','d','e')
        h.queue.DOWNLOAD_TARGET=2;h.queue.PREFETCH_LIMIT=4
        s=h.scheduler();s.step()
        self.assertEqual(len(s.downloads),2)
        self.assertEqual(s.prefetch_limit(),4)

    def test_실행중분리샘플은_선다운로드수에서제외(self):
        h=self.harness('orphan','ready1','ready2','next')
        for name in ['orphan','ready1','ready2']:h.success(name)
        h.queue.MAX_ACTIVE_SAMPLES=1
        h.queue.recovery_stage_running=lambda request: bool(request.get('wait_for_stage'))
        h.save(h.base/'작업재개요청.json',{'samples':{'orphan':{'stage':'10','requested_at':1,'wait_for_stage':{'pid':123}}}})
        s=h.scheduler();s.step()
        self.assertIn('next',s.downloads)
        self.assertNotIn('orphan',{r['sample_id'] for r in s.prefetch_rows()})
        self.assertEqual(len(s.prefetch_rows()),3)

    def test_완료즉시_같은cycle에_준비샘플투입(self):
        h=self.harness('first','second');h.queue.MAX_ACTIVE_SAMPLES=1
        h.success('first');h.success('second');s=h.scheduler();s.step()
        h.done.add('first');s.active['first'][1].result=0;s.step()
        self.assertIn('second',s.active)
        event=[detail for name,detail in h.events if name=='processing_ready_sample' and detail['sample']=='second'][-1]
        self.assertEqual(event['slot_released_to_start_seconds'],0)

    def test_고정슬롯없이_자원검사로_다섯째투입(self):
        names = [f'ready{i}' for i in range(6)]
        h = self.harness(*names)
        h.queue.MAX_ACTIVE_SAMPLES = None
        h.queue.resource_check = lambda active, candidate=None, prefetch=(): (
            len(active) < 5, 'RAM 예약 초과' if len(active) >= 5 else '')
        for name in names:
            h.success(name)
        s = h.scheduler()
        s.step()
        self.assertEqual(len(s.active), 5)
        queue_event = [detail for name, detail in h.events if name == 'sample_ready_queue'][-1]
        self.assertEqual(queue_event['ready_waiting'], ['ready5'])
        self.assertIsNone(queue_event['max_compute_slots'])
        self.assertIn('RAM 예약 초과', queue_event['waiting'][0]['reason'])

    def test_실행준비실패는_다른샘플막지않음(self):
        h=self.harness('bad','healthy');h.success('bad');h.success('healthy')
        def prepare(tag,*args):
            if tag=='1':raise RuntimeError('의도한 준비 오류')
        h.queue.prepare=prepare;s=h.scheduler();s.step()
        self.assertIn('bad',s.ledger['workers']);self.assertIn('healthy',s.active)

    def test_다운로드생성실패도_격리하고다음진행(self):
        h=self.harness('bad','healthy');start=h.queue.start
        def broken(kind,tag,*args,**kwargs):
            if kind=='download-worker' and tag=='1':raise OSError('의도한 생성 실패')
            return start(kind,tag,*args,**kwargs)
        h.queue.start=broken;s=h.scheduler();s.step()
        self.assertIn('healthy',s.downloads)
        self.assertEqual(s.ledger['downloads']['bad']['attempts'],1)
        self.assertGreater(s.ledger['downloads']['bad']['next_retry_at'],h.clock)

    def test_선다운로드3개는_소비즉시다음자료보충(self):
        names=['a','b','c','d','ready1','ready2','ready3','next']
        h=self.harness(*names)
        for name in names[:-1]:h.success(name)
        s=h.scheduler();s.step();self.assertIn('next',s.downloads)
        h.done.add('a');s.active['a'][1].result=0;s.step()
        self.assertIn('ready1',s.active);self.assertIn('next',s.downloads)
        self.assertEqual(len(s.prefetch_rows()),3)

    def test_계산4개면_다운로드4개를_선행실행(self):
        names=[f'ready{i}' for i in range(4)]+[f'new{i}' for i in range(6)]
        h=self.harness(*names)
        for name in names[:4]: h.success(name)
        s=h.scheduler();s.step()
        self.assertEqual(len(s.active),4)
        self.assertEqual(len(s.downloads),4)
        self.assertEqual(len(s.prefetch_rows()),4)
        event=[detail for name,detail in h.events if name=='sample_ready_queue'][-1]
        self.assertEqual(event['download_target'],4)
        self.assertEqual(event['prefetch_limit'],8)

    def test_준비버퍼가_차도_계산투입후_다운로드재개(self):
        names=[f'ready{i}' for i in range(12)]+['new']
        h=self.harness(*names)
        for name in names[:12]: h.success(name)
        s=h.scheduler();s.step()
        self.assertEqual(len(s.active),4)
        self.assertFalse(s.downloads)
        h.done.add('ready0');s.active['ready0'][1].result=0;s.step()
        self.assertIn('new',s.downloads)

    def test_제어주기1초(self):
        from unittest.mock import patch
        h=self.harness('x');s=h.scheduler()
        with patch.object(s,'step',side_effect=[None,0]),patch('handoff.ready_queue.time.sleep') as sleep:
            self.assertEqual(s.run(),0)
        sleep.assert_called_once_with(1.0)

    def test_분리단계중에는_중복재개안함_종료후재개(self):
        h = self.harness('orphan', 'normal')
        h.success('orphan'); h.success('normal')
        alive = {'value': True}
        h.queue.recovery_stage_running = lambda request: bool(request.get('wait_for_stage')) and alive['value']
        h.save(h.base / '작업재개요청.json', {'samples': {'orphan': {'stage': '10', 'requested_at': 1, 'wait_for_stage': {}}}})
        request = json.loads((h.base / '작업재개요청.json').read_text())
        request['samples']['orphan']['wait_for_stage'] = {'pid': 123}
        h.save(h.base / '작업재개요청.json', request)
        s = h.scheduler(); s.step()
        self.assertNotIn('orphan', s.active)
        self.assertIn('normal', s.active)
        self.assertEqual({r['sample_id'] for r in s.active_rows()}, {'orphan', 'normal'})
        alive['value'] = False
        s.step()
        self.assertIn('orphan', s.active)

    def test_분리단계도_동시샘플상한에포함(self):
        h = self.harness('orphan', 'normal')
        h.success('orphan'); h.success('normal')
        h.queue.MAX_ACTIVE_SAMPLES = 1
        h.queue.recovery_stage_running = lambda request: bool(request.get('wait_for_stage'))
        h.save(h.base / '작업재개요청.json', {'samples': {'orphan': {'stage': '10', 'requested_at': 1, 'wait_for_stage': {'pid': 123}}}})
        s = h.scheduler(); s.step()
        self.assertFalse(s.active)
        self.assertEqual([r['sample_id'] for r in s.active_rows()], ['orphan'])

    def test_다운로드실패중_준비샘플우선실행(self):
        h = self.harness('ready', 'failed', 'next')
        h.success('ready')
        h.records['failed'] = {'status': 'failed', 'finished': 0, 'error': 'HTTP 403'}
        s = h.scheduler()
        self.assertIsNone(s.step())
        self.assertEqual(h.started[0][:2], ('worker', 'ready'))
        self.assertIn('failed', s.downloads)
        s.downloads['failed'][1].result = 1
        s.step()
        self.assertIn('ready', s.active)
        self.assertGreater(s.ledger['downloads']['failed']['next_retry_at'], h.clock)
        self.assertIn('next', s.downloads)

    def test_재시도시간횟수_재시작보존(self):
        h = self.harness('retry')
        s = h.scheduler()
        s.step()
        s.downloads['retry'][1].result = 1
        s.step()
        count = s.ledger['downloads']['retry']['attempts']
        next_time = s.ledger['downloads']['retry']['next_retry_at']
        resumed = h.scheduler()
        resumed.step()
        self.assertFalse(resumed.downloads)
        self.assertEqual(resumed.ledger['downloads']['retry']['attempts'], count)
        h.clock = next_time + 1
        resumed.step()
        self.assertEqual(resumed.ledger['downloads']['retry']['attempts'], count + 1)

    def test_영구실패_첫세묶음밖_후속진행(self):
        h = self.harness('bad1', 'bad2', 'bad3', 'healthy')
        s = h.scheduler()
        for sample in ['bad1', 'bad2', 'bad3']:
            s.ledger['downloads'][sample] = {'attempts': 3, 'next_retry_at': 0, 'last_error': 'HTTP', 'deferred': False}
        s.step()
        self.assertIn('healthy', s.downloads)
        self.assertTrue(all(s.ledger['downloads'][x]['deferred'] for x in ['bad1', 'bad2', 'bad3']))

    def test_미완료남으면_완료보고안함(self):
        h = self.harness('bad', 'done')
        h.done.add('done')
        s = h.scheduler()
        s.ledger['downloads']['bad'] = {'attempts': 3, 'next_retry_at': 0, 'last_error': 'HTTP', 'deferred': True}
        self.assertEqual(s.step(), 1)
        self.assertEqual(h.events[-1][0], 'complete_with_pending')
        self.assertEqual(h.events[-1][1]['remaining_samples'], ['bad'])

    def test_완료샘플_원본없어도재다운로드안함(self):
        h = self.harness('done')
        h.done.add('done')
        self.assertEqual(h.scheduler().step(), 0)
        self.assertFalse(any(kind.startswith('download') for kind, _, _ in h.started))

    def test_자원부족_계산다운로드둘다대기(self):
        h = self.harness('ready', 'pending')
        h.success('ready')
        h.allowed = False
        s = h.scheduler()
        self.assertIsNone(s.step())
        self.assertFalse(h.started)

    def test_계산실패격리_정상계산계속(self):
        h = self.harness('bad', 'healthy')
        h.success('bad'); h.success('healthy')
        s = h.scheduler(); s.step()
        s.active['bad'][1].result = 1
        s.step()
        self.assertIn('healthy', s.active)
        self.assertIn('bad', s.ledger['workers'])

    def test_기존다운로드잠금_중복생성안함(self):
        h = self.harness('downloading')
        s = h.scheduler(); ctx = s.context(0)
        key = (str(h.queue.download_state(ctx['state'], ctx['wave'])), 'downloading')
        h.locks.add(key)
        s.adopt_workers(); s.step()
        self.assertIsInstance(s.downloads['downloading'][1], ExistingDownload)
        self.assertFalse(h.started)
        h.success('downloading'); h.locks.remove(key)
        s.step()
        self.assertIn('downloading', s.active)

    def test_다운로드완료후_다른묶음계산즉시투입(self):
        h = self.harness('slow', 'fast')
        s = h.scheduler(); s.step()
        h.success('fast'); s.downloads['fast'][1].result = 0
        s.step()
        self.assertIn('fast', s.active)
        self.assertIn('slow', s.downloads)

    def test_명시적복구요청_우선실행하고반복실패시보류(self):
        h = self.harness('healthy', 'recover')
        h.success('healthy'); h.success('recover')
        s = h.scheduler(); s.queue = h.queue
        h.queue.MAX_ACTIVE_SAMPLES = 1
        s.ledger['workers']['recover'] = {'reason': '이전 실패'}
        h.save(h.base / '작업재개요청.json', {'samples': {'recover': {'stage': '12', 'requested_at': 123}}})
        s.step()
        self.assertIn('recover', s.active)
        s.active['recover'][1].result = 1
        s.step()
        self.assertIn('recover', s.ledger['workers'])
        self.assertIn('healthy', s.active)



if __name__ == '__main__':
    unittest.main(verbosity=2)
