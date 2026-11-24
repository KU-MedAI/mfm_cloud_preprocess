"""같은 sample ID의 별도 CPU/GPU 실행을 출력 경로로 구분함."""
from pathlib import Path
from functools import lru_cache


@lru_cache(maxsize=4096)
def _read_environment(pid,started):
    try:
        raw=(Path('/proc')/str(pid)/'environ').read_bytes().split(b'\0')
        return {k:v for item in raw if b'=' in item for k,v in [item.decode().split('=',1)]
                if k in {'OUTPUT_ROOT','RUN_ID'}}
    except (OSError,UnicodeError):return {}


def _process_environment(pid):
    if not pid:return {}
    try:
        value=(Path('/proc')/str(int(pid))/'stat').read_text()
        started=int(value[value.rfind(')')+2:].split()[19])
    except (OSError,ValueError):return {}
    return _read_environment(int(pid),started)


def lease_matches_config(lease,config):
    """과거 장부에 소유 경로가 없으면 중복 계산 방지를 위해 보수적으로 포함함."""
    if not config or not config.get('output_root'):return True
    owner={}
    if lease.get('output_root'):
        owner={'OUTPUT_ROOT':lease['output_root'],'RUN_ID':lease.get('run_id')}
    else:
        for key in ['manager_pid','child_pid']:
            owner=_process_environment(lease.get(key))
            if owner.get('OUTPUT_ROOT'):break
    if not owner.get('OUTPUT_ROOT'):return True
    if Path(owner['OUTPUT_ROOT']).resolve()!=Path(config['output_root']).resolve():return False
    # batch_config가 이미 output_root에 실행 하위 경로를 반영함.
    # 실제 출력이 같으면 run_id가 달라도 같은 산출물을 보호해야 함.
    return True
