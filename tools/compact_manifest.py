"""새 산출물과 검증된 재사용 산출물의 완료 염색체를 공통 집계한다."""
from __future__ import annotations

import json
from pathlib import Path


def completed_contigs(manifest: dict, layout: Path, sample: str) -> set[str]:
    """재사용은 저장된 full 검증을 확인하며 독립 trace 검사는 호출자가 유지한다."""
    records = manifest.get("records", {})
    if not isinstance(records, dict):
        raise ValueError("compact records는 객체여야 함")
    completed: set[str] = set()
    for record in records.values():
        if not isinstance(record, dict):
            raise ValueError("compact record는 객체여야 함")
        if record.get("status") not in {"pass", "skip_existing"}:
            continue
        contig = record.get("contig")
        if (not isinstance(contig, str) or contig in {"", ".", ".."}
                or Path(contig).name != contig or contig in completed):
            raise ValueError(f"compact contig 누락/중복/잘못된 경로: {contig}")
        if record.get("status") == "skip_existing":
            directory = Path(layout) / contig
            if (record.get("sample") != sample
                    or not record.get("zarr")
                    or Path(record["zarr"]).resolve() != directory.resolve()):
                raise ValueError(f"재사용 compact sample/경로 불일치: {contig}")
            metadata = json.loads((directory / "contig_metadata.json").read_text())
            validation = metadata.get("validation", {})
            provenance = metadata.get("input_provenance", {})
            by_hap = validation.get("by_haplotype", {})
            fingerprint = provenance.get("fingerprint_sha256", "")
            full = (
                metadata.get("sample_id") == sample and metadata.get("contig") == contig
                and metadata.get("validation_level") == "full"
                and validation.get("level") == "full"
                and all(by_hap.get(hap, {}).get("status") == "pass"
                        and by_hap.get(hap, {}).get("level") == "full"
                        for hap in ("HAP1", "HAP2"))
                and provenance.get("contig") == contig
                and provenance.get("validation_level") == "full"
                and isinstance(fingerprint, str) and len(fingerprint) == 64
                and all(char in "0123456789abcdef" for char in fingerprint)
            )
            if not full:
                raise ValueError(f"재사용 compact의 full 검증/입력 provenance 불일치: {contig}")
        completed.add(contig)
    return completed
