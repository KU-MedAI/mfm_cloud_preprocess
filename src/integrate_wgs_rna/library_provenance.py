"""프로젝트 Methods에 근거한 WGS/RNA library 증폭 provenance.

플랫폼 이름이나 duplicate flag로 PCR 여부를 추론하지 않는다. 이 표는 현재 동결된
production cohort에만 적용하며, 새 프로젝트는 근거를 확인한 뒤 명시적으로 추가한다.
"""
from __future__ import annotations

from typing import Mapping


PROVENANCE_COLUMNS = (
    "wgs_library_amplification", "wgs_library_protocol", "wgs_library_evidence",
    "rna_library_amplification", "rna_library_protocol", "rna_library_pcr_cycles",
    "rna_library_evidence",
)
WGS_AMPLIFICATION_VALUES = {
    "pcr_free_confirmed", "pcr_free_protocol", "pcr_free_inferred", "pcr_amplified_confirmed",
}
RNA_AMPLIFICATION_VALUES = WGS_AMPLIFICATION_VALUES | {"not_applicable"}

WGS_PROJECT_PROVENANCE: dict[str, dict[str, str]] = {
    "1KGP_ILLUMINA": {
        "wgs_library_amplification": "pcr_free_confirmed",
        "wgs_library_protocol": "Illumina TruSeq DNA PCR-Free High Throughput Library Prep Kit",
        "wgs_library_evidence": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9439720/",
    },
    "HPRC_R1": {
        "wgs_library_amplification": "pcr_free_protocol",
        "wgs_library_protocol": "PacBio SMRTbell Express Template Prep Kit 2.0; size-selected HMW DNA",
        "wgs_library_evidence": "https://www.nature.com/articles/s41586-023-05896-x",
    },
    "HPRC_R2": {
        "wgs_library_amplification": "pcr_free_inferred",
        "wgs_library_protocol": "PacBio HiFi WGS standard library; sample-level PCR field unavailable",
        "wgs_library_evidence": "https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables/sequencing_data",
    },
    "HGSVC2": {
        "wgs_library_amplification": "pcr_free_protocol",
        "wgs_library_protocol": "PacBio SMRTbell/HiFi Express; sheared and size-selected HMW DNA",
        "wgs_library_evidence": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8026704/",
    },
    "HGSVC3": {
        "wgs_library_amplification": "pcr_free_protocol",
        "wgs_library_protocol": "PacBio SMRTbell Prep v3/Express 3.0 or prior HGSVC protocol",
        "wgs_library_evidence": "https://www.nature.com/articles/s41586-025-09140-6",
    },
}
GEUVADIS_RNA_PROVENANCE = {
    "rna_library_amplification": "pcr_amplified_confirmed",
    "rna_library_protocol": "Illumina TruSeq RNA Sample Prep Kit v2; adapter-ligated cDNA enrichment",
    "rna_library_pcr_cycles": "10",
    "rna_library_evidence": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1010958",
}


class LibraryProvenanceError(ValueError):
    """Library provenance가 허용된 계약과 다른 경우."""


def _clean(value: object) -> str:
    text = str(value or "").strip()
    return "" if text == "-" else text


def project_manifest_fields(project_id: str, has_rna: bool) -> dict[str, str]:
    """현재 production project의 canonical manifest 값을 반환한다."""
    if project_id not in WGS_PROJECT_PROVENANCE:
        raise LibraryProvenanceError(f"등록되지 않은 production project provenance: {project_id}")
    values = dict(WGS_PROJECT_PROVENANCE[project_id])
    values.update(GEUVADIS_RNA_PROVENANCE if has_rna else {
        "rna_library_amplification": "not_applicable", "rna_library_protocol": "",
        "rna_library_pcr_cycles": "", "rna_library_evidence": "",
    })
    return values


def validate_library_provenance(row: Mapping[str, object], *, has_rna: bool) -> dict[str, object]:
    """한 manifest 행의 provenance를 검사하고 JSON용 중첩 구조로 정규화한다."""
    project_id = _clean(row.get("project_id"))
    actual = {column: _clean(row.get(column)) for column in PROVENANCE_COLUMNS}
    if project_id in WGS_PROJECT_PROVENANCE:
        expected = project_manifest_fields(project_id, has_rna)
        mismatches = [f"{column}={actual[column]!r} expected={expected[column]!r}"
                      for column in PROVENANCE_COLUMNS if actual[column] != expected[column]]
        if mismatches:
            raise LibraryProvenanceError(
                f"project library provenance 불일치: project={project_id}; " + "; ".join(mismatches))
    elif not any(actual.values()):
        return {
            "wgs": {"amplification": "unspecified", "protocol": None, "evidence": None},
            "rna": {"amplification": "unspecified" if has_rna else "not_applicable",
                    "protocol": None, "pcr_cycles": None, "evidence": None},
        }
    else:
        if not all(actual[column] for column in (
                "wgs_library_amplification", "wgs_library_protocol", "wgs_library_evidence")):
            raise LibraryProvenanceError("사용자 project의 WGS provenance 세 필드는 모두 필요합니다")
        if has_rna and not all(actual[column] for column in (
                "rna_library_amplification", "rna_library_protocol", "rna_library_evidence")):
            raise LibraryProvenanceError("RNA pair의 RNA provenance와 evidence가 필요합니다")
        if not has_rna and actual["rna_library_amplification"] != "not_applicable":
            raise LibraryProvenanceError("RNA가 없는 sample은 rna_library_amplification=not_applicable이어야 합니다")

    if actual["wgs_library_amplification"] not in WGS_AMPLIFICATION_VALUES:
        raise LibraryProvenanceError(
            f"허용되지 않은 WGS amplification 상태: {actual['wgs_library_amplification']}")
    if actual["rna_library_amplification"] not in RNA_AMPLIFICATION_VALUES:
        raise LibraryProvenanceError(
            f"허용되지 않은 RNA amplification 상태: {actual['rna_library_amplification']}")
    cycles = actual["rna_library_pcr_cycles"]
    if actual["rna_library_amplification"] == "pcr_amplified_confirmed":
        if not cycles.isdigit() or int(cycles) <= 0:
            raise LibraryProvenanceError("PCR 증폭 RNA에는 양의 rna_library_pcr_cycles가 필요합니다")
    elif cycles:
        raise LibraryProvenanceError("PCR 증폭이 아닌 RNA에는 pcr_cycles를 기록할 수 없습니다")

    return {
        "wgs": {"amplification": actual["wgs_library_amplification"],
                "protocol": actual["wgs_library_protocol"], "evidence": actual["wgs_library_evidence"]},
        "rna": {"amplification": actual["rna_library_amplification"],
                "protocol": actual["rna_library_protocol"] or None,
                "pcr_cycles": int(cycles) if cycles else None,
                "evidence": actual["rna_library_evidence"] or None},
    }
