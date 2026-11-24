"""실제 quota와 구분한 정상·추가 논리 CPU의 공통 한도."""
import os
from .pool_leases import regular_cpu_claim, logical_cpu_claim
def loan_limits(total, setting):
    """실제 quota와 별개인 논리·대여 한도를 모든 입장 경로에서 같이 계산한다."""
    claims = max(2, min(4, float(setting.get('max_claim_pools', 2))))
    loans = max(1, min(3, float(setting.get('max_loan_pools', 1))))
    if 'max_extra_cpus' in setting:
        extra = max(0, min(total, int(setting['max_extra_cpus'])))
        return min(int(total * claims), total + extra), min(int(total * loans), extra)
    return int(total * claims), int(total * loans)




def available_loan_cpus(total, setting, leases, *, requested=None):
    """자기 lease를 제외한 예약에서 확보 가능한 총 대여 CPU임. 풀 잠금이 필요함."""
    if setting.get('enabled') is not True or os.environ.get('MFM_FIXED_SAMPLE_TOKEN'):
        return 0
    leases = tuple(leases)
    claims, loans = loan_limits(total, setting)
    claimed = sum(logical_cpu_claim(row) for row in leases)
    borrowed = sum(logical_cpu_claim(row) for row in leases if row.get('cpu_backfill'))
    capacity = max(0, min(claims - claimed, loans - borrowed))
    return capacity if requested is None else min(capacity, max(0, int(requested)))


def available_regular_cpus(total, setting, leases):
    """정상 CPU 예산과 전체 논리 상한을 함께 지킴. 호출자 자신의 lease는 제외함."""
    leases = tuple(leases)
    regular = sum(regular_cpu_claim(row) for row in leases)
    claimed = sum(logical_cpu_claim(row) for row in leases)
    ceiling, _ = loan_limits(total, setting)
    return max(0, min(total - regular, ceiling - claimed))
