"""WhatsHap 실행 전 donor alias와 전체 BAM/reference 범위를 검증한다.

원본 @RG는 보존한다. Manifest에 명시된 단일 donor alias만 read-group
무시를 허용하며, --regions로 read를 잘라내는 복구 방법은 사용하지 않는다.
"""
from __future__ import annotations

from pathlib import Path

import pysam


def validate_haplotag_contract(
    alignment: Path, tag_vcf: Path, reference: Path, sample: str,
    source_sample_id: str | None = None,
) -> dict[str, object]:
    alias = (source_sample_id or sample).strip()
    if not alias:
        raise ValueError("source_sample_id는 비어 있을 수 없습니다")
    with pysam.VariantFile(str(tag_vcf)) as variants:
        vcf_samples = list(variants.header.samples)
    if vcf_samples != [sample]:
        raise ValueError(f"단일 donor VCF가 필요합니다: expected={sample}, observed={vcf_samples}")
    with pysam.AlignmentFile(str(alignment), reference_filename=str(reference)) as bam:
        header = bam.header.to_dict()
        read_groups = header.get("RG", [])
        observed = {str(row.get("SM", "")).strip() for row in read_groups}
        if not read_groups or "" in observed:
            raise ValueError("BAM의 모든 read-group에 donor SM이 필요합니다")
        if len(observed) != 1 or not observed.issubset({sample, alias}):
            raise ValueError(f"단일 donor/명시적 alias만 허용됩니다: {sorted(observed)}; allowed={[sample, alias]}")
        # WhatsHap 기본 경로는 BAM 전체 contig와 unmapped read를 보존한다.
        # autosome-only FASTA를 넘기면 BAM chrX 처리 중 뒤늦게 실패할 수 있다.
        with pysam.FastaFile(str(reference)) as fasta:
            lengths = dict(zip(fasta.references, fasta.lengths))
        missing = [name for name in bam.references if name not in lengths]
        mismatch = [name for name, length in zip(bam.references, bam.lengths)
                    if name in lengths and lengths[name] != length]
        if missing or mismatch:
            raise ValueError(f"WhatsHap 전체 BAM reference 불일치: missing={missing}, length_mismatch={mismatch}")
    return {
        "sample": sample,
        "source_sample_id": alias,
        "bam_sm": sorted(observed),
        "vcf_samples": vcf_samples,
        "ignore_read_groups": sample not in observed,
        "reference_fasta": str(reference.resolve()),
        "bam_contig_count": len(header.get("SQ", [])),
        "read_preservation": "all_bam_contigs_and_unmapped; no_regions_restriction",
    }
