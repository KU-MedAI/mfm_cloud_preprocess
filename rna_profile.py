"""RNA 정렬·coverage profile 의 단일 정의 + 재사용 판정용 fingerprint.

15–19와 21이 이 값을 import해서 쓴다. 전체 fingerprint는 19의 수동 `check` 검증에 쓰고,
정렬 전용 fingerprint는 15·16 캐시에 써 coverage 정책만 바뀔 때 STAR index·정렬을 다시 만들지 않는다.

기준: GTEx/TOPMed V11 RNA-seq pipeline (broadinstitute/gtex-pipeline) — reference 가 HAP 개인 서열이고 annotation 이
HAP 좌표 GENCODE 인 것만 다르다.
"""
from __future__ import annotations

import hashlib
import os
import json
from pathlib import Path

PROFILE_NAME = "gtex_v11_hap"

# 16 — GTEx rnaseq/src/run_STAR.py 기본값 (26–57줄). 입출력·스레드·RG 는 16 이 붙인다.
STAR_GENOMIC_LENGTH_LOG2SCALE = -0.25        # STAR·GTEx 기본(GTEx 는 지정 안 함). 21이 HAP 비교 전에 이 항을 뺀다
STAR_ALIGN_ARGS = [
    "--twopassMode", "Basic", "--outSAMtype", "BAM", "SortedByCoordinate",
    "--outSAMattributes", "NH", "HI", "AS", "nM", "NM",
    "--outSAMunmapped", "Within", "--outSAMstrandField", "intronMotif",
    "--outFilterType", "BySJout", "--outFilterMultimapNmax", "20",
    "--outFilterMismatchNmax", "999", "--outFilterMismatchNoverLmax", "0.1",
    "--outFilterScoreMinOverLread", "0.33", "--outFilterMatchNminOverLread", "0.33",
    "--limitSjdbInsertNsj", "1200000", "--outFilterIntronMotifs", "None",
    "--alignIntronMin", "20", "--alignIntronMax", "1000000", "--alignMatesGapMax", "1000000",
    "--alignSJoverhangMin", "8", "--alignSJDBoverhangMin", "1",
    "--alignSoftClipAtReferenceEnds", "Yes", "--genomeLoad", "NoSharedMemory",
    "--quantMode", "GeneCounts", "--scoreGenomicLengthLog2scale", str(STAR_GENOMIC_LENGTH_LOG2SCALE),
]

# 15 — index 에 HAP chr1–22 와 함께 넣는 contig: GTEx 필터(TOPMed_RNAseq_pipeline.md 65–70: _alt·HLA·_decoy 제외) 후 chr1–22 를 뺀 것
from integrate_wgs_rna.processing_buffers import star_arguments
STAR_ALIGN_ARGS += star_arguments(os.environ)

INDEX_EXTRA_CONTIG_RULE = "Broad Homo_sapiens_assembly38 minus (*_alt, HLA*, *_decoy, chr1-22)"

# 17 — GTEx rnaseq/src/run_MarkDuplicates.py 39–48줄 (입출력·TMP_DIR·metrics 경로는 17 이 붙인다)
MARKDUP_ARGS = ["PROGRAM_RECORD_ID=null", "MAX_RECORDS_IN_RAM=500000", "SORTING_COLLECTION_SIZE_RATIO=0.25",
                "ASSUME_SORT_ORDER=coordinate", "TAGGING_POLICY=DontTag", "OPTICAL_DUPLICATE_PIXEL_DISTANCE=100"]

# 18 — GTEx bam2coverage.py의 `samtools view -q 255 -F 2816`와 같은 선택 정책이다.
# samtools depth 기본 필터에서 DUP를 해제하려면 `-g 1024`가 필요하다.
COVERAGE_DEPTH_ARGS = ["-a", "-J", "-Q", "255", "-G", "2820", "-g", "1024"]
COVERAGE_POLICY = "MAPQ 255, primary, QC pass, duplicate 포함(-g 1024), splice N 제외, deletion 포함, mate 겹침 2회"
COVERAGE_POLICY_REVISION = "2026-10-02_duplicate_included"
COVERAGE_DUPLICATE_POLICY = "included"


def coverage_metadata() -> dict:
    """실제 계산 인자와 같은 GTEx-compatible coverage 정책을 기록한다."""
    return {
        "read_filter": "GTEx bam2coverage equivalent: MAPQ255 primary QC-pass; duplicate included; supplementary excluded",
        "coverage_policy": COVERAGE_POLICY,
        "coverage_policy_revision": COVERAGE_POLICY_REVISION,
        "duplicate_policy": COVERAGE_DUPLICATE_POLICY,
        "effective_duplicate_policy": COVERAGE_DUPLICATE_POLICY,
        "legacy_filter_metadata_mismatch": False,
        "samtools_depth_args": list(COVERAGE_DEPTH_ARGS),
    }


def validate_coverage_metadata(attrs) -> None:
    """구형 duplicate-excluded coverage를 재사용하지 않는다."""
    expected = coverage_metadata()
    valid = all(attrs.get(k) == v for k, v in expected.items())
    if not valid:
        raise ValueError("coverage source 불일치: 실효 duplicate/filter policy")

# 14 — lift 정책
LIFT_POLICY = "trim"


def select_nonautosome(names: list[str]) -> list[str]:
    autosomes = {f"chr{i}" for i in range(1, 23)}
    return [n for n in names if n not in autosomes and not (n[-4:] == "_alt" or n[:3] == "HLA" or n[-6:] == "_decoy")]


def sha256_file(path: str | Path, chunk: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def fingerprint(gencode_gtf: str, extra_contig_source_fai: str = "", addnl: bool = True,
                align_only: bool = False) -> dict:
    """전체 또는 STAR 정렬 전용 profile fingerprint를 계산한다."""
    parts = {
        "profile": PROFILE_NAME,
        "star_align_args": STAR_ALIGN_ARGS,
        "lift_policy": LIFT_POLICY,
        "gencode_sha256": sha256_file(gencode_gtf),
        "index_extra_contigs": ({"rule": INDEX_EXTRA_CONTIG_RULE,
                                 "source_fai_sha256": sha256_file(extra_contig_source_fai) if extra_contig_source_fai else None}
                                if addnl else None),
    }
    if os.environ.get('MFM_RNA_BACKEND','cpu')=='parabricks':
        from integrate_wgs_rna.gpu_rna import effective_contract
        parts['gpu_alignment_contract']=effective_contract()
    if not align_only:
        parts["markdup_args"] = MARKDUP_ARGS
        parts["coverage_depth_args"] = COVERAGE_DEPTH_ARGS
    parts["sha256"] = hashlib.sha256(json.dumps(parts, sort_keys=True).encode()).hexdigest()
    return parts


if __name__ == "__main__":        # ... rna_profile.py <gencode_gtf> <source.fai|-> <0|1> [align]
    import sys
    g, fai, add = sys.argv[1], sys.argv[2], sys.argv[3] == "1"
    print(fingerprint(g, "" if fai == "-" else fai, add,
                      align_only=len(sys.argv) > 4 and sys.argv[4] == "align")["sha256"])
