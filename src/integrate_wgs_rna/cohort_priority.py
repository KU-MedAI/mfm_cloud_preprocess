"""WGS-RNA pair 우선 cohort 선정기.

로컬 후보 TSV의 WGS/RNA 파일과 인덱스 검증 결과를 다시 확인하고,
WGS-RNA pair donor를 WGS-only보다 먼저 고정한다.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import random
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

PAIR = "wgs_rna_pair"
WGS_ONLY = "wgs_only"
SELECTION_ORDER = (PAIR, WGS_ONLY)
PLATFORMS = ("ont", "illumina", "hifi")
REQUIRED_COLUMNS = (
    "sample_id",
    "dataset_role",
    "platform",
    "alignment_path",
    "alignment_exists",
    "alignment_has_index",
    "rna_path",
    "rna_exists",
    "rna_has_index",
)


class CohortPriorityError(ValueError):
    """후보 TSV 또는 선정 설정이 계약을 위반했을 때 발생한다."""


def _present(value: str | None) -> bool:
    return bool(value and value.strip() and value.strip() != "-")


def _truthy(value: str | None) -> bool:
    return (value or "").strip().lower() in {"1", "true", "yes", "y"}


def normalise_platform(value: str | None) -> str | None:
    key = (value or "").strip().lower().replace("-", "_").replace(" ", "")
    aliases = {
        "ont": "ont",
        "nanopore": "ont",
        "oxfordnanopore": "ont",
        "oxford_nanopore": "ont",
        "illumina": "illumina",
        "illumina_30x": "illumina",
        "ilmn": "illumina",
        "hifi": "hifi",
        "hifi_ccs": "hifi",
        "pacbiohifi": "hifi",
        "pacbio_hifi": "hifi",
        "pacbioccs": "hifi",
        "pacbio_ccs": "hifi",
        "ccs": "hifi",
    }
    return aliases.get(key)


def selection_role(row: dict[str, str]) -> str | None:
    role = (row.get("dataset_role") or "").strip().lower()
    use = (row.get("recommended_training_use") or "").strip().lower()
    if role == "post_training" or use == "post_training_wgs_rna_pair":
        return PAIR
    if role == "pre_training" or use == "pre_training_wgs_only":
        return WGS_ONLY
    return None


def _index_candidates(path: Path) -> tuple[Path, ...]:
    suffix = path.suffix.lower()
    candidates = [Path(f"{path}.crai"), Path(f"{path}.bai")]
    if suffix == ".cram":
        candidates.append(path.with_suffix(".crai"))
    elif suffix == ".bam":
        candidates.append(path.with_suffix(".bai"))
    return tuple(dict.fromkeys(candidates))


def _file_failure(
    row: dict[str, str],
    prefix: str,
    *,
    require_index: bool,
    check_paths: bool,
) -> str | None:
    path_text = row.get(f"{prefix}_path")
    if not _present(path_text):
        return f"{prefix}_path_missing"
    if not _truthy(row.get(f"{prefix}_exists")):
        return f"{prefix}_file_not_verified"
    if require_index and not _truthy(row.get(f"{prefix}_has_index")):
        return f"{prefix}_index_not_verified"

    path = Path(str(path_text))
    if check_paths and not path.is_file():
        return f"{prefix}_file_missing"
    if require_index and check_paths and not any(index.is_file() for index in _index_candidates(path)):
        return f"{prefix}_index_missing"
    return None


def _candidate_failure(
    row: dict[str, str],
    role: str,
    *,
    check_paths: bool,
) -> str | None:
    platform = normalise_platform(row.get("platform"))
    if platform is None:
        return "unsupported_platform"
    failure = _file_failure(row, "alignment", require_index=True, check_paths=check_paths)
    if failure:
        return failure
    if role == PAIR:
        return _file_failure(row, "rna", require_index=True, check_paths=check_paths)
    return None


def default_targets() -> dict[str, dict[str, int]]:
    return {
        PAIR: {"ont": 0, "illumina": 5, "hifi": 5},
        WGS_ONLY: {"ont": 0, "illumina": 5, "hifi": 5},
    }


def select_candidates(
    rows: list[dict[str, str]],
    targets: dict[str, dict[str, int]],
    *,
    seed: int = 42,
    check_paths: bool = True,
) -> dict[str, Any]:
    """pair-first·donor-unique 정책으로 후보를 선택한다."""

    pools: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    rejections: list[dict[str, str]] = []

    for row in rows:
        role = selection_role(row)
        platform = normalise_platform(row.get("platform"))
        if role is None:
            reason = "unsupported_dataset_role"
        elif platform is None:
            reason = "unsupported_platform"
        else:
            reason = _candidate_failure(row, role, check_paths=check_paths)
        if reason:
            rejections.append({
                "sample_id": row.get("sample_id", ""),
                "selection_role": role or "-",
                "normalized_platform": platform or "-",
                "reason": reason,
                "alignment_path": row.get("alignment_path", ""),
                "rna_path": row.get("rna_path", ""),
            })
            continue
        assert role is not None and platform is not None
        pools[(role, platform)].append(row)

    rng = random.Random(seed)
    selected: list[dict[str, str]] = []
    used_samples: set[str] = set()
    summary: list[dict[str, Any]] = []

    for role in SELECTION_ORDER:
        for platform in PLATFORMS:
            target = int(targets.get(role, {}).get(platform, 0))
            if target < 0:
                raise CohortPriorityError(f"target은 0 이상이어야 합니다: {role}/{platform}={target}")
            pool = sorted(
                pools.get((role, platform), []),
                key=lambda row: (
                    row.get("sample_id", ""),
                    row.get("alignment_path", ""),
                    row.get("rna_path", ""),
                ),
            )
            rng.shuffle(pool)
            chosen_count = 0
            for row in pool:
                sample = row.get("sample_id", "").strip()
                if sample in used_samples:
                    rejections.append({
                        "sample_id": sample,
                        "selection_role": role,
                        "normalized_platform": platform,
                        "reason": "donor_already_selected",
                        "alignment_path": row.get("alignment_path", ""),
                        "rna_path": row.get("rna_path", ""),
                    })
                    continue
                if chosen_count >= target:
                    rejections.append({
                        "sample_id": sample,
                        "selection_role": role,
                        "normalized_platform": platform,
                        "reason": "quota_reached",
                        "alignment_path": row.get("alignment_path", ""),
                        "rna_path": row.get("rna_path", ""),
                    })
                    continue
                accepted = dict(row)
                accepted["selection_role"] = role
                accepted["normalized_platform"] = platform
                accepted["selection_order"] = str(len(selected) + 1)
                selected.append(accepted)
                used_samples.add(sample)
                chosen_count += 1

            summary.append({
                "selection_role": role,
                "platform": platform,
                "target": target,
                "eligible": len(pool),
                "selected": chosen_count,
                "missing": target - chosen_count,
            })

    rejection_counts = dict(sorted(Counter(item["reason"] for item in rejections).items()))
    qc = {
        "schema": "mfm_cohort_priority_v1",
        "seed": seed,
        "check_paths": check_paths,
        "selection_order": list(SELECTION_ORDER),
        "targets": targets,
        "summary": summary,
        "selected_total": len(selected),
        "selected_unique_donors": len(used_samples),
        "rejection_counts": rejection_counts,
    }
    return {"selected": selected, "rejections": rejections, "qc": qc}


def read_candidate_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise CohortPriorityError(f"TSV header가 없습니다: {path}")
        headers = list(reader.fieldnames)
        missing = [column for column in REQUIRED_COLUMNS if column not in headers]
        if missing:
            raise CohortPriorityError("필수 후보 열이 없습니다: " + ",".join(missing))
        rows = list(reader)
    return headers, rows


def _atomic_csv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(fd, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({column: row.get(column, "") for column in columns})
        os.replace(temporary, path)
    except Exception:
        Path(temporary).unlink(missing_ok=True)
        raise


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        Path(temporary).unlink(missing_ok=True)
        raise


def write_selection(
    output_dir: Path,
    headers: list[str],
    result: dict[str, Any],
) -> dict[str, Path]:
    """선정 결과를 pair/WGS-only manifest와 QC로 원자적으로 기록한다."""

    extra = ["selection_role", "normalized_platform", "selection_order"]
    selected_columns = headers + [column for column in extra if column not in headers]
    rejection_columns = [
        "sample_id",
        "selection_role",
        "normalized_platform",
        "reason",
        "alignment_path",
        "rna_path",
    ]
    paths = {
        "selected": output_dir / "selected.tsv",
        "pair_manifest": output_dir / "pair_manifest.tsv",
        "wgs_only_manifest": output_dir / "wgs_only_manifest.tsv",
        "rejections": output_dir / "selection_rejections.tsv",
        "qc": output_dir / "selection_qc.json",
    }
    selected = result["selected"]
    _atomic_csv(paths["selected"], selected_columns, selected)
    _atomic_csv(
        paths["pair_manifest"],
        selected_columns,
        [row for row in selected if row["selection_role"] == PAIR],
    )
    _atomic_csv(
        paths["wgs_only_manifest"],
        selected_columns,
        [row for row in selected if row["selection_role"] == WGS_ONLY],
    )
    _atomic_csv(paths["rejections"], rejection_columns, result["rejections"])
    _atomic_json(paths["qc"], result["qc"])
    return paths


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="WGS-RNA pair를 우선하는 local cohort selector")
    parser.add_argument("--candidates", type=Path, required=True, help="local 파일 상태를 포함한 후보 TSV")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--pair-ont", type=int, default=0)
    parser.add_argument("--pair-illumina", type=int, default=5)
    parser.add_argument("--pair-hifi", type=int, default=5)
    parser.add_argument("--wgs-ont", type=int, default=0)
    parser.add_argument("--wgs-illumina", type=int, default=5)
    parser.add_argument("--wgs-hifi", type=int, default=5)
    parser.add_argument(
        "--no-path-check",
        action="store_true",
        help="이미 검증된 snapshot TSV만 재현할 때 실제 파일·index 재확인을 생략한다",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    headers, rows = read_candidate_tsv(args.candidates)
    targets = default_targets()
    targets[PAIR].update({"ont": args.pair_ont, "illumina": args.pair_illumina, "hifi": args.pair_hifi})
    targets[WGS_ONLY].update({"ont": args.wgs_ont, "illumina": args.wgs_illumina, "hifi": args.wgs_hifi})
    result = select_candidates(rows, targets, seed=args.seed, check_paths=not args.no_path_check)
    paths = write_selection(args.output_dir, headers, result)
    print(f"selected_total\t{result['qc']['selected_total']}")
    for item in result["qc"]["summary"]:
        print(
            f"{item['selection_role']}\t{item['platform']}\t"
            f"selected={item['selected']}\ttarget={item['target']}\tmissing={item['missing']}"
        )
    print(f"pair_manifest\t{paths['pair_manifest']}")
    print(f"wgs_only_manifest\t{paths['wgs_only_manifest']}")
    print(f"selection_qc\t{paths['qc']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
