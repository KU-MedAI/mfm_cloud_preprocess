#!/usr/bin/env python3
"""HG03125 active source를 PacBio CLR에서 동일 개인 1KGP Illumina 30x로 교체한다."""
from __future__ import annotations

import csv
import json
import os
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SAMPLE = "HG03125"
CRAM_URL = "https://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989126/HG03125.final.cram"
CRAI_URL = CRAM_URL + ".crai"
CRAM_BYTES = "15575174940"


def read(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    tmp = path.with_name(f".{path.name}.part.{os.getpid()}")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)
    os.replace(tmp, path)


def update_preprocess() -> None:
    path = ROOT / "manifests/preprocess_manifest.tsv"
    fields, rows = read(path)
    template = next(r for r in rows if r["sample_id"] == "HG00100")
    target = next(r for r in rows if r["sample_id"] == SAMPLE)
    keep = {
        "sample_id": SAMPLE, "dataset_role": target["dataset_role"], "data_group": "wgs",
        "rna_path": "", "rna_url": "", "rna_index_url": "", "rna_input_bytes": "0",
        "rna_input_gib": "0.000", "rna_index_bytes": "0",
        "rna_library_amplification": "not_applicable", "rna_library_protocol": "",
        "rna_library_pcr_cycles": "", "rna_library_evidence": "",
    }
    target.update(template)
    target.update(keep)
    target.update({
        "platform": "illumina", "read_type": "short_read",
        "project_id": "1KGP_ILLUMINA", "cohort": "1KGP_ILLUMINA",
        "alignment_path": f"wgs/1KGP_ILLUMINA/{SAMPLE}/{SAMPLE}.final.cram",
        "read_files_manifest": "", "wgs_url": CRAM_URL, "wgs_index_url": CRAI_URL,
        "wgs_input_bytes": CRAM_BYTES, "wgs_input_gib": "14.506",
        "known_input_total_bytes": CRAM_BYTES, "known_input_total_gib": "14.506",
        "size_completeness": "primary_only", "size_basis": "source_manifest_bytes",
        "public_vcf_sample_id": "", "public_vcf_path": "", "vcf_source": "illumina_panel",
    })
    write(path, fields, rows)

    rf = ROOT / "manifests/read_files.tsv"
    rf_fields, rf_rows = read(rf)
    write(rf, rf_fields, [r for r in rf_rows if r["sample_id"] != SAMPLE])


def update_catalog() -> None:
    path = ROOT / "share_package/samples_3212.tsv"
    fields, rows = read(path)
    template = next(r for r in rows if r["sample_id"] == "HG00100")
    target = next(r for r in rows if r["sample_id"] == SAMPLE)
    memberships = target.get("all_published_pacbio_memberships", "HGSVC2")
    role = target.get("recommended_training_use", "pre_training_wgs_only")
    target.update(template)
    target.update({
        "sample_id": SAMPLE, "selected_platform": "Illumina_30x",
        "selected_project": "1KGP_ILLUMINA",
        "all_published_pacbio_memberships": memberships,
        "read_type": "Illumina PCR-free 30x", "read_format": "GRCh38 aligned CRAM+CRAI",
        "read_alignment_state": "공개 aligned CRAM 사용 가능",
        "read_primary_url_or_rule": CRAM_URL,
        "wgs_rna_pair_status": "WGS_only", "has_rna_bam_strict": "False",
        "rna_project": "", "rna_file_format": "", "rna_reference": "",
        "rna_primary_url": "", "rna_index_url": "", "rna_primary_bytes": "0",
        "rna_manifest_status": "not_applicable",
        "wgs_primary_url": CRAM_URL, "wgs_index_url": CRAI_URL,
        "wgs_primary_bytes": CRAM_BYTES, "wgs_size_basis": "source_manifest_bytes",
        "recommended_training_use": role,
        "selection_reason": "CLR 제외 정책에 따라 동일 개인의 1KGP Illumina 30x를 선택",
        "package_version": "mfm_share_package_3212_20260929_v3",
    })
    write(path, fields, rows)
    write(ROOT / "share_package/illumina_2918_samples.tsv", fields,
          [r for r in rows if r["selected_platform"] == "Illumina_30x"])
    write(ROOT / "share_package/pacbio_294_samples.tsv", fields,
          [r for r in rows if r["selected_platform"] == "PacBio_HiFi"])
    write(ROOT / "share_package/rna_pair_449.tsv", fields,
          [r for r in rows if r["has_rna_bam_strict"] == "True"])

    csv_path = ROOT / "share_package/samples_3212.csv"
    if csv_path.exists():
        tmp = csv_path.with_name(f".{csv_path.name}.part.{os.getpid()}")
        with tmp.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
            writer.writeheader(); writer.writerows(rows)
        os.replace(tmp, csv_path)


def update_active_downloads() -> None:
    package = ROOT / "share_package"
    illumina = package / "manifests/illumina_30x_files.tsv"
    fields, rows = read(illumina)
    template = next(r for r in rows if r["sample_id"] == "HG00100")
    rows = [r for r in rows if r["sample_id"] != SAMPLE]
    target = dict(template)
    target.update({
        "sample_id": SAMPLE, "project": "1KGP_ILLUMINA", "file_url": CRAM_URL,
        "index_url": CRAI_URL, "size_bytes": CRAM_BYTES,
        "status": "url_identified_checksum_not_frozen",
    })
    rows.append(target)
    rows.sort(key=lambda r: r["sample_id"])
    write(illumina, fields, rows)

    direct = package / "manifests/wgs_direct_download.tsv"
    dfields, drows = read(direct)
    templates = [dict(r) for r in drows if r["sample_id"] == "HG00100" and r["dataset"] == "Illumina_30x"]
    drows = [r for r in drows if r["sample_id"] != SAMPLE]
    for row in templates:
        role = row["file_role"]
        row.update({
            "sample_id": SAMPLE, "project": "1KGP_ILLUMINA",
            "direct_url": CRAM_URL if role == "primary" else CRAI_URL,
            "expected_bytes": CRAM_BYTES if role == "primary" else "",
            "source_manifest": "illumina_2918_samples.tsv",
            "notes": "HG03125 CLR 제외; 동일 개인 1KGP Illumina 30x 선택",
        })
        drows.append(row)
    drows.sort(key=lambda r: (r["dataset"], r["sample_id"], r["file_role"], r["direct_url"]))
    write(direct, dfields, drows)


def update_phase_and_summary() -> None:
    package = ROOT / "share_package"
    path = package / "phase_assignments.tsv"
    fields, rows = read(path)
    target = next(r for r in rows if r["sample_id"] == SAMPLE)
    target.update({
        "estimated_download_bytes": CRAM_BYTES,
        "estimated_download_gib": f"{int(CRAM_BYTES) / 1024**3:.3f}",
        "selected_platform": "Illumina_30x", "selected_project": "1KGP_ILLUMINA",
        "wgs_download_status": "url_identified_checksum_not_frozen",
        "estimated_wgs_bytes": CRAM_BYTES, "estimated_rna_bytes": "0",
        "wgs_size_basis": "source_manifest_bytes",
    })
    write(path, fields, rows)

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["phase"]].append(row)
    phase_fields = [
        "phase", "samples", "estimated_download_bytes", "estimated_download_tib",
        "strict_rna_bam_pair", "PacBio_HiFi", "PacBio_CLR", "Illumina_30x",
        "contains_estimated_pacbio_sizes",
    ]
    phase_rows = []
    for phase, members in sorted(grouped.items()):
        total = sum(int(r["estimated_download_bytes"] or 0) for r in members)
        platforms = Counter(r["selected_platform"] for r in members)
        phase_rows.append({
            "phase": phase, "samples": str(len(members)), "estimated_download_bytes": str(total),
            "estimated_download_tib": f"{total / 1024**4:.3f}",
            "strict_rna_bam_pair": str(sum(r["has_rna_bam_strict"] == "True" for r in members)),
            "PacBio_HiFi": str(platforms["PacBio_HiFi"]), "PacBio_CLR": str(platforms["PacBio_CLR"]),
            "Illumina_30x": str(platforms["Illumina_30x"]),
            "contains_estimated_pacbio_sizes": str(any("estimated" in r["wgs_size_basis"] for r in members)),
        })
    write(package / "phase_summary.tsv", phase_fields, phase_rows)

    registry_fields, registry = read(package / "samples_3212.tsv")
    del registry_fields
    platforms = Counter(r["selected_platform"] for r in registry)
    projects = Counter(r["selected_project"] for r in registry)
    counts = [
        ("cohort", "unique_samples", len(registry)), ("cohort", "in_1kgp_3202_panel", 3202),
        ("cohort", "outside_1kgp_panel", 10),
        ("training", "strict_rna_bam_pair", sum(r["has_rna_bam_strict"] == "True" for r in registry)),
        ("training", "wgs_only", sum(r["has_rna_bam_strict"] != "True" for r in registry)),
    ]
    counts.extend(("platform", key, value) for key, value in sorted(platforms.items()))
    counts.extend(("project", key, value) for key, value in sorted(projects.items()))
    pair_platform = Counter(r["selected_platform"] for r in registry if r["has_rna_bam_strict"] == "True")
    counts.extend(("rna_pair_platform", key, value) for key, value in sorted(pair_platform.items()))
    cpath = package / "summary_counts.tsv"
    tmp = cpath.with_name(f".{cpath.name}.part.{os.getpid()}")
    with tmp.open("w") as handle:
        handle.write("category\tlabel\tcount\n")
        for category, label, count in counts:
            handle.write(f"{category}\t{label}\t{count}\n")
    os.replace(tmp, cpath)

    pacbio_ids = {r["sample_id"] for r in registry if r["selected_platform"] == "PacBio_HiFi"}
    _, pacbio_files = read(package / "manifests/pacbio_file_manifest.tsv")
    active_read_files = sum(r["sample_id"] in pacbio_ids for r in pacbio_files)
    summary = {
        "package_schema": "mfm_share_package_3212_20260929_v3", "generated_date": "2026-09-29",
        "cohort": {
            "unique_samples": 3212, "in_1kgp_3202_panel": 3202, "outside_1kgp_panel": 10,
            "platform": dict(sorted(platforms.items())), "project": dict(sorted(projects.items())),
        },
        "rna": {
            "strict_bam_pairs": 449, "wgs_only": 2763,
            "pair_platform": dict(sorted(pair_platform.items())),
        },
        "readiness": {
            "illumina_file_urls_identified": platforms["Illumina_30x"],
            "pacbio_donors_frozen_active": platforms["PacBio_HiFi"],
            "pacbio_read_files_frozen_active": active_read_files,
            "historical_pacbio_freeze_donors": 295,
            "historical_pacbio_read_files": len(pacbio_files),
            "public_vcf_files_frozen": 704, "reference_files_frozen": 3,
        },
        "historical_predecessor": "../260826/3202_sample_share_package",
    }
    spath = package / "summary.json"
    tmp = spath.with_name(f".{spath.name}.part.{os.getpid()}")
    tmp.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, spath)


def main() -> int:
    update_preprocess()
    update_catalog()
    update_active_downloads()
    update_phase_and_summary()
    print("HG03125: PacBio_CLR -> 1KGP Illumina_30x; active sample count unchanged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
