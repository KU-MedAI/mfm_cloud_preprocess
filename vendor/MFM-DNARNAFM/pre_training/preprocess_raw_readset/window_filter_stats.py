"""창 유지 조건을 바꾸지 않고 blacklist 우선으로 제외 사유를 집계한다."""
from __future__ import annotations

import json
from decimal import Decimal
from pathlib import Path
from typing import Any

DROP_REASON_PRIORITY = "blacklist_then_n"
N_FRACTION_DENOMINATOR = "total_windows_before_filters"


def exclusion_reason(*, ok_n: bool, ok_blacklist: bool) -> str:
    if not ok_blacklist:
        return "blacklist"
    if not ok_n:
        return "n"
    return ""


def index_filter_policy(path: Path, rows: list[dict[str, Any]]) -> dict[str, Any]:
    """구형 TSV의 중첩률은 반드시 당시 실행의 임계값과 비교한다."""
    meta_path = path.with_name("windows.json")
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    needs_legacy_policy = any(
        "drop_reason" not in row
        and str(row.get("kept", "")).strip().lower() not in {"1", "true"}
        and int(float(row.get("n_nonacgt", "0") or 0)) > 0
        for row in rows
    )
    if needs_legacy_policy and "max_blacklist_frac" not in meta:
        raise ValueError(f"구형 N 포함 제외 창 재집계에는 당시 max_blacklist_frac가 있는 windows.json이 필요합니다: {path}")
    return {
        "max_blacklist_frac": float(meta.get("max_blacklist_frac", .01)),
        "drop_n": bool(meta.get("drop_n", True)),
    }


def count_index_exclusions(
    rows: list[dict[str, Any]], *, max_blacklist_frac: float, drop_n: bool = True,
) -> dict[str, int]:
    counts = {"dropped_blacklist": 0, "dropped_n": 0}
    threshold = Decimal(str(max_blacklist_frac))
    for row in rows:
        kept = str(row.get("kept", "")).strip().lower() in {"1", "true"}
        if "drop_reason" in row:
            reason = str(row["drop_reason"]).strip()
            if reason not in {"", "blacklist", "n"} or kept != (reason == ""):
                raise ValueError("windows.tsv의 kept와 drop_reason이 불일치합니다")
        elif kept:
            reason = ""
        elif not drop_n or int(float(row.get("n_nonacgt", "0") or 0)) == 0:
            reason = "blacklist"
        else:
            if row.get("bl_frac", "") == "":
                raise ValueError("구형 N 포함 제외 창의 blacklist 우선 재집계에 bl_frac가 필요합니다")
            fraction = Decimal(str(row["bl_frac"]))
            # 과거 TSV는 소수 6자리 반올림값이다. 경계에서는 원 판정을 추정하지 않는다.
            half_step = Decimal("0.0000005")
            if fraction.as_tuple().exponent == -6 and fraction-half_step <= threshold <= fraction+half_step:
                raise ValueError("구형 bl_frac가 반올림 경계입니다. 원 FASTA/chain으로 window index를 재검증해 drop_reason을 기록해야 합니다")
            reason = "blacklist" if fraction > threshold else "n"
        if reason:
            counts["dropped_" + reason] += 1
    return counts
