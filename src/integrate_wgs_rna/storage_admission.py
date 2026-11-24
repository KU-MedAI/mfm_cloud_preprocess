"""현재 디스크 여유로 대형 입력 준비 동시성을 자동 산정한다."""
from pathlib import Path
import os
import time


def refresh_storage_capacity(policy,total):
    setting=policy.get('adaptive_prepare_storage',{})
    if setting.get('enabled') is not True:
        return None
    folder=setting.get('path',policy.get('storage_path'))
    per_job=int(setting.get('peak_bytes_per_job',350_000_000_000))
    reserve=int(setting.get('reserve_bytes',150_000_000_000))
    if not folder or per_job<=0 or reserve<0:
        raise ValueError('입력 준비 디스크 정책이 올바르지 않음')
    stat=os.statvfs(folder);free=stat.f_bavail*stat.f_frsize
    slots=max(0,min(total//4,max(0,free-reserve)//per_job))
    previous=policy.setdefault('stage_concurrency',{}).get('WGS:01')
    policy['stage_concurrency']['WGS:01']=slots
    changed=previous!=slots or 'prepare_storage_plan' not in policy
    report={'schema':'mfm_adaptive_prepare_storage_v1','time':time.time(),
            'path':str(folder),'free_bytes':free,'reserve_bytes':reserve,
            'peak_bytes_per_job':per_job,'previous_slots':previous,'slots':slots,'policy_changed':changed}
    # 매주기 timestamp만으로 정책을 다시 쓰지 않는다.
    if changed:
        policy['prepare_storage_plan']=report
    return report
