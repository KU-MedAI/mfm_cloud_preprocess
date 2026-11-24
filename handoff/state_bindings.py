"""이전 실행의 데이터 장부만 연결함. 외부 운영 Python은 실행·import하지 않음."""
from pathlib import Path
from contextlib import contextmanager
import copy
import fcntl
import hashlib
import json


def read(path):
    p = Path(path)
    return json.loads(p.read_text()) if p.is_file() else {}


class StateBindings:
    def __init__(self, config, rows):
        self.config = config
        self.rows = {r['sample_id']: r for r in rows}
        self.path = Path(config['state_root']) / 'imported_samples.json'
        data = read(self.path)
        self.source = Path(data['source']) if data.get('source') else None
        if data and data.get('schema') != 'mfm_imported_sample_states_v1':
            raise ValueError('이관 장부 schema가 다릅니다.')
        self.items = data.get('samples', {})
        self.mirrored = {}
        if set(self.items) - set(self.rows):
            raise ValueError('이관 장부와 현재 선택 범위가 다릅니다.')
        for sample, item in self.items.items():
            path = Path(item['selection'])
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != item['selection_sha256']:
                raise ValueError('이관된 선택 파일이 변경됨: ' + sample)
            old = read(item['config'])
            for key in ['data_root', 'output_root']:
                if Path(old[key]).resolve() != Path(config[key]).resolve():
                    raise ValueError('이관 작업의 입력/출력 경로가 다름: ' + sample)
            if Path(old['state_root']).resolve() != Path(item['state_root']).resolve():
                raise ValueError('이관 작업 상태 경로가 다름: ' + sample)

    @contextmanager
    def control_lock(self):
        """이전 운영 진입점을 다시 실행해도 동시 제어되지 않게 같은 잠금을 사용함."""
        path = self.source / '연속운영.lock' if self.source else None
        if path is None or not path.exists():
            yield
            return
        with path.open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            yield

    def state(self, sample):
        return Path(self.items.get(sample, {}).get('state_root', self.config['state_root']))

    def download_state(self, sample, default):
        return Path(self.items.get(sample, {}).get('download_state_root', default))

    def record(self, sample):
        return read(self.state(sample) / 'samples' / (sample + '.json'))

    def launch(self, sample, kind, recovery, api, directory):
        """이전 task index·작업 경로를 유지하고 새 공유 자원/복구 정책만 전달함."""
        item = self.items[sample]
        config = copy.deepcopy(read(item['config']))
        # 새 제어 큐의 정책을 사용하되 선택 fingerprint/run_id/state는 유지함.
        config.update(cpu_pool_limit=self.config.get('cpu_pool_limit', 'auto'),
                      memory_limit_bytes=self.config['memory_limit_bytes'])
        config['cpu_allocation_mode']=self.config.get('cpu_allocation_mode','dynamic')
        for key in ['sample_cpus', 'local_jobs', 'resource_resolution', 'resources_resolved']:
            config[key] = copy.deepcopy(self.config[key])
        config.setdefault('pipeline_env', {}).update(self.config.get('pipeline_env', {}))
        request = recovery or {}
        env = config.setdefault('sample_pipeline_env', {}).setdefault(sample, {})
        env.update(self.config.get('sample_pipeline_env', {}).get(sample, {}))
        if request:
            env['WGS_START_AT'] = request['stage']
            if request.get('overwrite_prepared'):
                env['STAGE04_OVERWRITE'] = '1'
            if request.get('rebuild_haplotag_partial'):
                env['HAPLOTAG_REBUILD_PARTIAL'] = '1'
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        if kind == 'download-worker':
            from .download_only import configure, initialize
            from .planner import read_tsv
            config = configure(config)
            if str(Path(config['state_root']).resolve()) != str(Path(item['download_state_root']).resolve()):
                raise ValueError('이관 다운로드 식별자가 달라짐: ' + sample)
            if not (Path(config['state_root']) / 'effective_config.json').is_file():
                initialize(config, read_tsv(item['selection']))
        elif not (Path(config['state_root']) / 'references.json').is_file():
            from .planner import read_tsv
            api.prepare(config, read_tsv(item['selection']))
        path = directory / (sample + '.' + kind + '.json')
        api.atomic_json(path, config)
        return path, Path(item['selection']), int(item['task_index'])

    def mirror(self, api):
        """통합 집계용 복사본만 만듦. 활성 worker는 이전 장부에 계속 기록함."""
        directory = Path(self.config['state_root']) / 'samples'
        for sample in self.items:
            source = self.state(sample) / 'samples' / (sample + '.json')
            if not source.is_file():
                continue
            stat = source.stat()
            fingerprint = (stat.st_mtime_ns, stat.st_size, stat.st_ino)
            if self.mirrored.get(sample) == fingerprint:
                continue
            record = self.record(sample)
            if record:
                api.atomic_json(directory / (sample + '.json'), record)
                self.mirrored[sample] = fingerprint
