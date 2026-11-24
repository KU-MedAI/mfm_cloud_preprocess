#!/usr/bin/env python3
"""27 — v2 통합 전처리 레이아웃용 sample catalog 생성.

    python 27_make_catalog.py <final_cohort_root> [--out PATH] [--meta-tsv PATH]

예:
    python 27_make_catalog.py \
      /path/to/output/wgs_rna_pair/1KGP_ILLUMINA

★ 왜 새로 만드나 (기존 `process_code/91_make_catalog.py` 를 안 쓰고)
  91 은 **구 파이프라인 전용**이다. 파일명 규약(`<S>.ont.hg38.haplotagged.bam`,
  `<S>.ont.personal.HAP1.cram`)과 RNA 경로(`raw/<S>/rna_1bp/total.zarr`)가 v2 와 다르고,
  코호트 단위 `sequence_variant_policy` 를 config.sh 에서 하나로 읽는다. 그런데 같은 코호트
  디렉터리 안에 **구 112명(snp_only)과 신규 10명(snv_indel)이 섞여 있다.** 코호트 단위 정책으로는
  이 둘을 구분할 수 없다.

  그래서 이 스크립트는 **샘플별 `metadata/sample.json` 을 단일 출처**로 삼는다.
  전처리가 스스로 적은 값이므로 추측이 없고, 구·신 혼재도 자연스럽게 구분된다.
  91 은 구 코호트용으로 그대로 남겨 둔다.

  배경: `post-training/journal/2609/post_training_data_issue_260903.md` §7

★ 원칙
  경로는 **존재를 확인하고** 적는다. 없는 것을 있다고 적으면 카탈로그가 있는 게 더 위험하다.
  `qc_all_paths_exist` 는 "이 샘플로 학습이 되는가"를 뜻하고, 판정 기준은 정책마다 다르다
  (snv_indel 신규 경로는 compact 이 곧 모델 입력이라 compact 존재가 핵심이다).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from tools.compact_manifest import completed_contigs
from catalog_common import (
    first_existing, load_population_meta, load_splits, read_clair3_qc, read_json,
 )


COLS = [
    "sample_id", "cohort", "dataset_role", "read_type", "project_id",
    "sequence_variant_policy", "platform", "reference_build", "sample_status",
    "layout_name", "hap_read_mode", "split",
    # 집단·가족 — split 을 짤 때 가족 누수를 막으려면 반드시 봐야 한다
    "sex", "pop", "superpop", "family_id", "father", "mother", "related",
    # 원본
    "source_wgs", "source_vcf",
    # 개인 서열
    "individual_fasta_hap1", "individual_fasta_hap2",
    "n_fasta_hap1", "n_fasta_hap2", "n_chain",
    # 정렬
    "haplotagged_bam", "taglist",
    "realigned_hap1", "realigned_hap2", "realign_target_hap1", "realign_target_hap2",
    # window index
    "windows_tsv", "windows_blacklist_tsv",
    # ★ 모델 입력 (v2 신규)
    "compact_root", "compact_status", "compact_contigs", "compact_validation_level",
    "trace_status", "trace_contigs",
    # ★ RNA (v2 신규 — tracks/ 가 아니라 rna_personal 아래에 있다)
    "rna_personal_root", "phase_mode",
    "coverage_zarr_hap1", "coverage_zarr_hap2",
    "ase_training_ready",
    "allelic_counts", "allelic_counts_unit",
    # phasing 품질 — phased_het_fraction 만으로는 파편화를 못 잡는다
    "phased_het_fraction", "phase_set_count", "phase_block_n50_bp",
    # 요약 플래그
    "has_rna", "has_compact", "has_realigned", "has_chain",
    "qc_all_paths_exist", "processing_status", "training_eligible", "missing",
    "final_validation_status", "final_validation_failed_checks",
]

CHROMS = [f"chr{i}" for i in range(1, 23)]

# haplotagged BAM 은 파이프라인 세대에 따라 이름이 다르다. 둘 다 받아준다.
HAPLOTAG_PATTERNS = (
    "{s}.clair3_longread.haplotagged.bam",
    "{s}.public_vcf_longread.haplotagged.bam",
    "{s}.ont.hg38.haplotagged.bam",
)
TAGLIST_PATTERNS = ("{s}.clair3_longread.taglist.tsv.gz", "{s}.taglist.tsv.gz")
REALIGN_PATTERNS = ("{s}.personal.{H}.cram", "{s}.ont.personal.{H}.cram",
                    "{s}.personal.{H}.bam", "{s}.ont.personal.{H}.bam")


def scan_sample(sample_root: Path, cohort: str, pop_meta: dict, split_of: dict) -> dict:
    s = sample_root.name
    meta = read_json(sample_root / "metadata" / "sample.json")
    hap_read_mode = meta.get("hap_read_mode", "all_to_both")
    policy = meta.get("sequence_variant_policy", "")
    role = meta.get("dataset_role", "")
    missing: list[str] = []

    h1, h2 = sample_root / "individual_fasta" / "HAP1", sample_root / "individual_fasta" / "HAP2"
    n1 = len(list(h1.glob("chr*.fa"))) if h1.is_dir() else 0
    n2 = len(list(h2.glob("chr*.fa"))) if h2.is_dir() else 0
    n_chain = (len(list(h1.glob("chr*.chain"))) + len(list(h2.glob("chr*.chain")))) if h1.is_dir() else 0

    hap_dir = sample_root / "haplotagged"
    bam = first_existing(hap_dir, HAPLOTAG_PATTERNS, s=s) if hap_dir.is_dir() else None
    tag = first_existing(hap_dir, TAGLIST_PATTERNS, s=s) if hap_dir.is_dir() else None

    realn: dict[str, str] = {}
    target: dict[str, str] = {}
    for hap in ("HAP1", "HAP2"):
        rdir = sample_root / "realigned" / hap_read_mode / hap
        if not rdir.is_dir():
            # 2026-09-27 이전 layout 호환. 신규 산출물은 mode 디렉터리가 필수다.
            rdir = sample_root / "realigned" / hap
        if not rdir.is_dir():
            continue
        found = first_existing(rdir, REALIGN_PATTERNS, s=s, H=hap)
        # 인덱스까지 있어야 완성으로 본다 — 절반짜리 경로를 적으면 카탈로그가 거짓말을 한다
        if found is None:
            continue
        idx = Path(str(found) + (".crai" if found.suffix == ".cram" else ".bai"))
        if not idx.exists():
            continue
        realn[hap] = str(found)
        tfa = rdir / "target.fa"
        target[hap] = str(tfa) if tfa.exists() else ""

    wtsv = sample_root / "windows" / "windows.tsv"
    wbl = sample_root / "windows" / "windows_blacklist.tsv"

    layout = meta.get("layout_name", "")
    compact_root = sample_root / "compact" / layout if layout else None
    compact_status, compact_contigs, compact_validation_level = "", 0, ""
    compact_pass_contigs: set[str] = set()
    if compact_root is not None and compact_root.is_dir():
        cman = read_json(compact_root / "manifest.json")
        compact_status = cman.get("status", "")
        compact_validation_level = cman.get("validation_level", "")
        try:
            compact_pass_contigs = completed_contigs(cman, compact_root, s)
        except (OSError, ValueError, TypeError, AttributeError):
            compact_status = "invalid_manifest"
            compact_pass_contigs = set()
        compact_contigs = len(compact_pass_contigs)
    else:
        compact_root = None

    trace_manifest = sample_root / "qc" / "stage12_trace" / hap_read_mode / "trace_manifest.json"
    trace = read_json(trace_manifest)
    trace_status = trace.get("status", "")
    trace_pass_contigs = {k for k, v in trace.get("contigs", {}).items() if v == "pass"}
    trace_contigs = len(trace_pass_contigs)
    trace_root_matches = bool(compact_root) and Path(trace.get("compact_root", "/")).resolve() == compact_root.resolve()

    rna = meta.get("rna_personal") or {}
    cov = rna.get("coverage_zarr", {}) if isinstance(rna, dict) else {}
    allelic = rna.get("allelic_counts", "") if isinstance(rna, dict) else ""
    unit = ""
    if allelic:
        side = Path(os.path.splitext(allelic)[0] + ".meta.json")
        unit = read_json(side).get("counting_unit", "") if side.is_file() else "read(legacy)"
    ase_qc_path = Path(rna.get("ase_qc_json", "")) if isinstance(rna, dict) and rna.get("ase_qc_json") else None
    ase_qc = read_json(ase_qc_path) if ase_qc_path else {}
    sample_qc = ase_qc.get("sample_qc", {})
    ase_training_ready = isinstance(sample_qc, dict) and sample_qc.get("sample_training_ready") is True

    qc = read_clair3_qc(sample_root)
    final_validation = read_json(sample_root / "qc" / "final_package_validation.json")
    final_validation_status = str(final_validation.get("status", ""))
    raw_failed_checks = final_validation.get("failed_checks")
    final_failed_checks = (
        (";".join(str(item) for item in raw_failed_checks) or "failed_checks_nonempty")
        if isinstance(raw_failed_checks, list) and raw_failed_checks else
        "" if isinstance(raw_failed_checks, list) else "failed_checks_schema_invalid"
    )
    final_pass = (final_validation_status == "pass"
                  and isinstance(raw_failed_checks, list) and not raw_failed_checks)
    pm = pop_meta.get(s, {})
    expected_contigs = set(CHROMS)
    has_compact = (
        bool(compact_root) and compact_status == "complete"
        and compact_validation_level == "full" and compact_pass_contigs == expected_contigs
    )
    has_trace = (
        trace_status == "pass" and trace_pass_contigs == expected_contigs and trace_root_matches
    )
    has_rna = (
        set(cov) == {"HAP1", "HAP2"}
        and all(Path(cov[hap]).is_dir() for hap in ("HAP1", "HAP2"))
        and bool(allelic) and Path(allelic).is_file()
        and ase_qc_path is not None and ase_qc_path.is_file()
        and bool(rna.get("ase_training_eligibility", ""))
        and Path(rna.get("ase_training_eligibility", "")).is_file()
    )
    has_realigned = len(realn) == 2
    has_chain = n_chain == 44

    if policy == "snv_indel":
        need = {"compact_full_22": has_compact, "source_trace_22": has_trace}
        if role == "post_training":
            need["rna_outputs"] = has_rna
    else:
        need = {"individual_fasta": n1 == 22 and n2 == 22,
                "haplotagged": bam is not None, "windows": wtsv.exists()}
    missing = [key for key, ok in need.items() if not ok]
    outputs_complete = not missing
    training_eligible = outputs_complete and final_pass and (role != "post_training" or ase_training_ready)
    if not final_pass:
        processing_status = "validation_pending" if not final_validation_status else "validation_failed"
    else:
        processing_status = ("complete" if training_eligible else
                             "processed_qc_fail" if outputs_complete else "incomplete")

    return {
        "sample_id": s, "cohort": cohort,
        "dataset_role": role, "read_type": meta.get("read_type", ""),
        "project_id": meta.get("project_id", cohort),
        "sequence_variant_policy": policy,
        "platform": meta.get("platform", ""), "reference_build": meta.get("reference_build", ""),
        "sample_status": meta.get("sample_status", ""),
        "layout_name": layout, "hap_read_mode": hap_read_mode,
        "split": split_of.get(s, ""),
        "sex": pm.get("sex", ""), "pop": pm.get("pop", ""), "superpop": pm.get("superpop", ""),
        "family_id": pm.get("family_id", ""), "father": pm.get("father", ""),
        "mother": pm.get("mother", ""), "related": pm.get("related", ""),
        "source_wgs": "", "source_vcf": "",          # source_manifest.tsv 에서 채운다(아래)
        "individual_fasta_hap1": str(h1) if n1 else "",
        "individual_fasta_hap2": str(h2) if n2 else "",
        "n_fasta_hap1": n1, "n_fasta_hap2": n2, "n_chain": n_chain,
        "haplotagged_bam": str(bam) if bam else "",
        "taglist": str(tag) if tag else "",
        "realigned_hap1": realn.get("HAP1", ""), "realigned_hap2": realn.get("HAP2", ""),
        "realign_target_hap1": target.get("HAP1", ""), "realign_target_hap2": target.get("HAP2", ""),
        "windows_tsv": str(wtsv) if wtsv.exists() else "",
        "windows_blacklist_tsv": str(wbl) if wbl.exists() else "",
        "compact_root": str(compact_root) if compact_root else "",
        "compact_validation_level": compact_validation_level,
        "trace_status": trace_status, "trace_contigs": trace_contigs,
        "compact_status": compact_status, "compact_contigs": compact_contigs,
        "rna_personal_root": rna.get("rna_personal_root", "") if isinstance(rna, dict) else "",
        "ase_training_ready": str(ase_training_ready),
        "phase_mode": rna.get("phase_mode", "") if isinstance(rna, dict) else "",
        "coverage_zarr_hap1": cov.get("HAP1", ""), "coverage_zarr_hap2": cov.get("HAP2", ""),
        "allelic_counts": allelic, "allelic_counts_unit": unit,
        "phased_het_fraction": qc.get("vcf_phased_het_fraction", ""),
        "phase_set_count": qc.get("vcf_phase_set_count", ""),
        "phase_block_n50_bp": qc.get("vcf_phase_block_n50_bp", ""),
        "has_rna": str(has_rna), "has_compact": str(has_compact),
        "has_realigned": str(has_realigned), "has_chain": str(has_chain),
        "qc_all_paths_exist": str(outputs_complete),
        "processing_status": processing_status,
        "training_eligible": str(training_eligible),
        "missing": ";".join(missing),
        "final_validation_status": final_validation_status,
        "final_validation_failed_checks": final_failed_checks,
    }


def fill_sources(row: dict, sample_root: Path) -> None:
    """metadata/source_manifest.tsv 에서 원본 경로를 가져온다."""
    path = sample_root / "metadata" / "source_manifest.tsv"
    if not path.is_file():
        return
    wgs, vcf = [], []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            kind, p = r.get("kind", ""), r.get("path", "")
            if not p:
                continue
            if kind == "source_wgs":
                wgs.append(p)
            elif kind == "source_vcf":
                vcf.append(p)
    row["source_wgs"] = ";".join(wgs)
    row["source_vcf"] = ";".join(vcf)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("final_cohort_root", help="final/<role>/<read_type>/<cohort>")
    ap.add_argument("--out", default="")
    ap.add_argument("--meta-tsv", default="",
                    help="선택 집단/가족 메타 TSV(sample, sex, pop, superpop, family_id, father, mother, related)")
    ap.add_argument("--splits-dir", default="",
                    help="split 정의 디렉터리의 부모(<...>/metadata). 기본값은 preprocessed_data/metadata")
    a = ap.parse_args()

    root = Path(a.final_cohort_root).resolve()
    if not root.is_dir():
        print(f"[FATAL] 없는 경로: {root}", file=sys.stderr)
        return 1
    cohort = root.name

    # split metadata는 기존 final layout이면 외부 metadata를 읽을 수 있지만,
    # catalog 자체는 언제나 전달 받은 cohort root 아래에 materialize한다.
    if a.splits_dir:
        pp_metadata = Path(a.splits_dir)
    elif len(root.parents) >= 4:
        pp_metadata = root.parents[3] / "metadata"
    else:
        # portable/handoff smoke root처럼 final hierarchy보다 얕은 경우에는
        # split metadata가 선택 사항이므로 root 아래를 비어 있는 기본값으로 쓴다.
        pp_metadata = root / "metadata"
    split_of = load_splits(pp_metadata)
    pop_meta = load_population_meta(Path(a.meta_tsv) if a.meta_tsv else None)

    # canonical final/<role>/<type>/<cohort>/<sample>뿐 아니라 handoff smoke
    # test의 <root>/wgs/<project>/<sample>도 지원한다. sample.json 하나당 한
    # sample만 취하고 catalog/ 같은 자체 산출물은 제외한다.
    sample_dirs = sorted({path.parent.parent for path in root.rglob("metadata/sample.json")
                          if "catalog" not in path.parts})
    if not sample_dirs:
        print(f"[FATAL] metadata/sample.json 을 가진 샘플이 없다: {root}", file=sys.stderr)
        return 1

    rows = []
    for d in sample_dirs:
        row = scan_sample(d, cohort, pop_meta, split_of)
        fill_sources(row, d)
        rows.append(row)

    out = Path(a.out) if a.out else root / "catalog" / "sample_catalog.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".part")
    with open(tmp, "w", newline="") as fh:
        # lineterminator 를 명시한다. csv 기본값 "\r\n" 이면 awk/split("\n") 쪽에서 마지막 컬럼에 \r 이 붙는다.
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t",
                           lineterminator="\n", extrasaction="raise")
        passed_rows = [row for row in rows if row["final_validation_status"] == "pass"
                       and not row["final_validation_failed_checks"]
                       and row["qc_all_paths_exist"] == "True"]
        w.writeheader(); w.writerows(passed_rows)
    os.replace(tmp, out)

    incomplete = [row for row in rows if row not in passed_rows]
    incomplete_path = out.parent / "incomplete_samples.tsv"
    tmp_incomplete = incomplete_path.with_suffix(incomplete_path.suffix + ".part")
    with open(tmp_incomplete, "w", newline="") as fh:
        fields = ("sample_id", "project_id", "dataset_role", "final_validation_status",
                  "failed_checks", "processing_status", "missing")
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in incomplete:
            writer.writerow({
                "sample_id": row["sample_id"], "project_id": row["project_id"],
                "dataset_role": row["dataset_role"],
                "final_validation_status": row["final_validation_status"],
                "failed_checks": row["final_validation_failed_checks"],
                "processing_status": row["processing_status"], "missing": row["missing"],
            })
    os.replace(tmp_incomplete, incomplete_path)

    print(f"[ok] {out}  (final pass {len(passed_rows)}/{len(rows)}행)")
    keys = ("sequence_variant_policy", "hap_read_mode", "compact_validation_level",
                "trace_status", "processing_status", "training_eligible",
                "has_compact", "has_rna", "has_realigned", "has_chain",
                "allelic_counts_unit", "qc_all_paths_exist")
    for label, selected_rows in (("발견한 전체 샘플", rows), ("최종 catalog 출력 샘플", passed_rows)):
        print(f"  [{label}] {len(selected_rows)}명")
        for key in keys:
            counts: dict[str, int] = {}
            for r in selected_rows:
                counts[str(r[key])] = counts.get(str(r[key]), 0) + 1
            print(f"    {key:24s} {dict(sorted(counts.items()))}")

    # 가족 누수 경고 — split 을 짜기 전에 반드시 봐야 한다
    related = [r["sample_id"] for r in rows if r.get("related") and r["related"] != "unrelated"]
    if related:
        print(f"  ⚠ 가족 관계 있는 샘플 {len(related)}명: {', '.join(related)}")
        for r in rows:
            if r["sample_id"] in related:
                print(f"      {r['sample_id']}: {r['related']} "
                      f"(family_id={r['family_id']} father={r['father']} mother={r['mother']})")

    bad = [r for r in rows if r["qc_all_paths_exist"] != "True"]
    if bad:
        print(f"[WARN] 필수 경로가 빠진 샘플 {len(bad)}/{len(rows)}")
        for r in bad[:20]:
            print(f"      {r['sample_id']}: missing={r['missing']}")
    else:
        print("  ✅ 전 샘플 필수 경로 존재")
    if incomplete:
        print(f"[WARN] final package validation/필수 출력 미충족 {len(incomplete)}명: {incomplete_path}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
