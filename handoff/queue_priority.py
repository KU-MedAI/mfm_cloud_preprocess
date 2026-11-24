"""기존 배치/작업 인덱스를 유지하며 read 종류별 선택 순서만 결정함."""
from pathlib import Path
import json


class PriorityPolicy:
    def __init__(self, scheduler):
        self.scheduler = scheduler
        self.spec = {}
        self.default_max_active = getattr(scheduler.q, 'MAX_ACTIVE_SAMPLES', None)
        self.refresh()

    def refresh(self):
        q = self.scheduler.q
        path = Path(q.BASE) / '우선처리정책.json'
        spec = json.loads(path.read_text()) if path.exists() else getattr(q, 'PRIORITY_POLICY', {})
        if not isinstance(spec, dict) or spec.get('mode', 'ready_first') not in {'ready_first', 'long_read_first', 'long_read_only'}:
            raise ValueError('우선처리 정책 mode는 ready_first, long_read_first 또는 long_read_only여야 함')
        drain = spec.get('drain_samples', [])
        blocked = spec.get('blocked_samples', {})
        if not isinstance(drain, list) or any(not isinstance(x, str) for x in drain):
            raise ValueError('drain_samples는 샘플ID 목록이어야 함')
        if not isinstance(blocked, dict) or any(not isinstance(k, str) or not isinstance(v, str) or not v for k, v in blocked.items()):
            raise ValueError('blocked_samples는 샘플ID:보류사유 매핑이어야 함')
        if not isinstance(spec.get('wait_for_drain', True), bool):
            raise ValueError('wait_for_drain은 bool이어야 함')
        if not isinstance(spec.get('pair_first', True), bool):
            raise ValueError('pair_first는 bool이어야 함')
        limit = spec.get('max_active_samples', self.default_max_active)
        if limit is not None and (isinstance(limit, bool) or not isinstance(limit, int) or limit < 1):
            raise ValueError('max_active_samples는 양의 정수 또는 null이어야 함')
        config=getattr(q,'config',{})
        if config.get('cpu_allocation_mode')=='fixed_sample':
            fixed_limit=int(config['local_jobs'])
            limit=fixed_limit if limit is None else min(limit,fixed_limit)
        allowed = spec.get('download_allowed_block_reasons', [])
        if not isinstance(allowed, list) or any(not isinstance(x, str) or not x for x in allowed):
            raise ValueError('download_allowed_block_reasons는 명시 보류사유 목록이어야 함')
        self.spec = spec
        # 실행 중 계산을 중단하지 않고 신규 입장부터 실시간 동시 수를 적용함.
        q.MAX_ACTIVE_SAMPLES = limit

    @property
    def enabled(self):
        return self.spec.get('mode') in {'long_read_first', 'long_read_only'}

    def paused(self, sample):
        if self.spec.get('mode') != 'long_read_only' or sample in self.spec.get('drain_samples', []):
            return False
        row = next((r for r in self.scheduler.rows if r['sample_id'] == sample), None)
        return row is not None and row.get('read_type') != 'long_read'

    def blocked(self, sample, *, worker=True):
        reason = self.spec.get('blocked_samples', {}).get(sample, '')
        if not worker and reason in self.spec.get('download_allowed_block_reasons', []):
            return ''
        return reason

    def pending_long(self):
        s = self.scheduler
        return [e for e in s.entries if e[2].get('read_type') == 'long_read'
                and e[2]['sample_id'] not in s.completed and not s.deferred(e)]

    def draining(self):
        s = self.scheduler
        return [e[2]['sample_id'] for e in s.entries
                if e[2]['sample_id'] in self.spec.get('drain_samples', [])
                and e[2]['sample_id'] not in s.completed
                and (e[2]['sample_id'] in s.active or not s.deferred(e))]

    def reason(self, row, *, worker):
        if not self.enabled:
            return ''
        if worker and self.spec.get('wait_for_drain', True) and self.draining():
            return '정책 전환 전 실행 샘플 완료 대기: ' + ','.join(self.draining())
        if row.get('read_type') != 'long_read' and (self.spec.get('mode') == 'long_read_only' or self.pending_long()):
            return 'long read 우선 처리 중: short read 신규 투입 대기'
        return ''

    def ordered(self, entries, *, worker):
        s = self.scheduler
        def recovery(e):
            sample = e[2]['sample_id']
            return (sample not in s.recoveries, s.recoveries.get(sample, {}).get('priority', 10))
        if not self.enabled:
            return sorted(entries, key=recovery) if worker else entries
        return sorted(entries, key=lambda e: (
            e[2].get('read_type') != 'long_read',
            e[2].get('data_group') != 'wgs_rna_pair' if e[2].get('read_type') == 'long_read' and self.spec.get('pair_first', True) else False,
            *recovery(e),
            int(e[2].get('known_input_total_bytes') or 0) if e[2].get('read_type') == 'long_read' else 0))

    def prefetch_count(self, candidate=None):
        s = self.scheduler
        rows = s.prefetch_rows()
        if self.enabled and (self.spec.get('mode') == 'long_read_only' or self.pending_long()):
            # 이미 준비된 short 파일은 보존하며 long 준비 버퍼를 차지하지 않음.
            rows = [r for r in rows if r.get('read_type') == 'long_read']
        if self.enabled and candidate and self.spec.get('pair_first', True) and candidate.get('data_group') == 'wgs_rna_pair':
            # 이미 확보한 WGS-only 원본 때문에 미준비 pair 입력을 막지 않는다.
            rows = [r for r in rows if r.get('data_group') == 'wgs_rna_pair']
        return len(rows)

    def status(self):
        return {'mode': self.spec.get('mode', 'ready_first'),
                'pair_first': self.spec.get('pair_first', True),
                'max_active_samples': self.scheduler.q.MAX_ACTIVE_SAMPLES,
                'drain_waiting': self.draining() if self.enabled and self.spec.get('wait_for_drain', True) else [],
                'preserved_active_samples': self.draining(),
                'wait_for_drain': self.spec.get('wait_for_drain', True),
                'long_remaining': len(self.pending_long()) if self.enabled else None,
                'blocked_samples': self.spec.get('blocked_samples', {}),
                'registered_long_reads': sum(r.get('read_type') == 'long_read' for r in self.scheduler.rows),
                'paused_short_samples': [r['sample_id'] for r in self.scheduler.rows
                    if r['sample_id'] not in self.scheduler.completed and self.paused(r['sample_id'])]}

