"""Stage 21 판정과 Stage 22의 소비 계약. 평균 rd는 설명용이지 통과 기준이 아니다."""
from __future__ import annotations

import math

CNV_POLICY_REVISION = "2026-10-01_per_hap_complete"


def evaluate_hap_depths(values: list[float | None], minimum: float, maximum: float) -> tuple[bool, str]:
    """모든 입력 HAP의 관측값이 있어야 판정한다. 결측은 CNV 정상의 근거가 아니다."""
    if not (math.isfinite(minimum) and math.isfinite(maximum) and 0 <= minimum <= maximum):
        raise ValueError("CNV rd 범위는 유한한 0 <= min <= max여야 합니다")
    if any(v is not None and (not math.isfinite(v) or v < 0) for v in values):
        raise ValueError("CNV rd는 유한한 음이 아닌 값이어야 합니다")
    if not values or any(v is None for v in values):
        return False, "missing_depth"
    passed = all(minimum <= v <= maximum for v in values)
    return passed, "pass" if passed else "outside_range"


def consume_cnv_row(row: dict[str, str] | None, *, enabled: bool) -> tuple[float | None, bool, str]:
    """평균으로 재판정하지 않고 producer pass_cnv를 따른다. 결측/비활성을 구분한다."""
    if not enabled:
        return None, True, "disabled"
    if row is None:
        return None, False, "missing_gene"
    if row.get("pass_cnv") not in {"0", "1"}:
        raise ValueError("CNV 입력에는 유효한 pass_cnv=0|1 판정이 필요합니다")
    rd = float(row["rd"]) if row.get("rd") not in (None, "", "NA") else None
    raw_haps = [row.get("rd_hap1"), row.get("rd_hap2")]
    observed = [float(v) for v in raw_haps if v not in (None, "", "NA")]
    if any(not math.isfinite(v) or v < 0 for v in observed + ([rd] if rd is not None else [])):
        raise ValueError("CNV 입력에 유효하지 않은 rd가 있습니다")
    if rd is None or len(observed) != 2:
        return rd, False, "missing_depth"
    passed = row["pass_cnv"] == "1"
    # producer 판정을 쓰므로 TSV 소수점 반올림으로 경계 판정이 뒤집히지 않는다.
    return rd, passed, "pass" if passed else "outside_range"
