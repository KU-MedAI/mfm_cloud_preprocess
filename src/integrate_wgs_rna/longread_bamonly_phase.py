from __future__ import annotations

import gzip
import json
import math
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, NamedTuple

from .cohort_vcf import expand_chroms

BASES = frozenset("ACGT")


@dataclass
class BamOnlyPhaseConfig:
    """BAM-only long-read read clustering configuration.

    v2 policy deliberately avoids Illumina/statistical VCF-based haplotagging for
    ONT/HiFi reads.  This module therefore discovers heterozygous-like SNV
    markers directly from the long-read BAM/CRAM, clusters reads by co-occurring
    allele patterns, and writes a QC/taglist contract.

    The cluster labels produced here are *not* HAP1/HAP2 labels yet.  Cluster
    orientation must be resolved later against HAP1/HAP2 personal FASTA, normally
    using personal-alignment score statistics.  Keeping this distinction explicit
    prevents downstream compact generation from silently treating arbitrary
    cluster numbers as biological HP tags.
    """

    sample_id: str
    alignment: Path
    output_dir: Path
    reference_fasta: Path | None = None
    contigs: str = "chr1-chr22"
    platform: str = "ont"
    min_mapq: int = 20
    min_base_quality: int = 0
    min_depth: int = 6
    min_marker_reads: int = 2
    min_minor_allele_fraction: float = 0.20
    max_major_allele_fraction: float = 0.80
    max_third_allele_fraction: float = 0.10
    max_marker_gap: int = 1_000_000
    max_markers_per_block: int = 5_000
    min_read_markers: int = 2
    min_score_delta: int = 2
    max_iterations: int = 6
    include_supplementary: bool = False
    exclude_duplicates: bool = False
    max_pileup_depth: int = 100_000
    marker_detection_chunk_bp: int = 1_000_000
    max_markers_per_contig: int = 0
    regions: str = ""
    unknown_read_policy: str = "include_both"
    dry_run: bool = False


class Region(NamedTuple):
    contig: str
    start: int | None
    end: int | None


@dataclass(frozen=True)
class Marker:
    contig: str
    pos0: int
    allele1: str
    allele2: str
    allele1_count: int
    allele2_count: int
    depth: int
    third_allele_fraction: float


@dataclass
class ClusterCall:
    read_name: str
    contig: str
    ps: int
    cluster: int
    marker_count: int
    informative_count: int
    score_cluster1: int
    score_cluster2: int
    score_delta: int
    confidence: float
    status: str


def run_bamonly_phase(config: BamOnlyPhaseConfig) -> dict[str, Any]:
    """Infer BAM-only read clusters and write marker/taglist/QC outputs.

    Output labels are intentionally named ``cluster`` rather than ``HP`` because
    this step does not yet orient clusters to HAP1/HAP2 personal FASTA.  The
    resulting taglist is suitable for QC and for a later orientation step; it is
    not a production replacement for a verified HP/PS-tagged alignment.
    """
    validate_config(config)
    contract = build_contract(config)
    if config.dry_run:
        return contract

    import pysam

    config.output_dir.mkdir(parents=True, exist_ok=True)
    marker_path = config.output_dir / f"{config.sample_id}.bamonly_phase.markers.tsv.gz"
    read_tag_path = config.output_dir / f"{config.sample_id}.bamonly_phase.read_clusters.tsv.gz"
    qc_path = config.output_dir / f"{config.sample_id}.bamonly_phase.qc.json"

    regions = parse_regions(config.regions) if config.regions else [Region(contig, None, None) for contig in expand_chroms([config.contigs])]
    qc: dict[str, Any] = {
        **contract,
        "status": "completed_unoriented_clusters",
        "outputs": {
            "markers_tsv_gz": str(marker_path),
            "read_clusters_tsv_gz": str(read_tag_path),
            "qc_json": str(qc_path),
        },
        "summary": {
            "regions_requested": len(regions),
            "regions_processed": 0,
            "markers_total": 0,
            "blocks_total": 0,
            "blocks_with_assignments": 0,
            "read_block_rows_total": 0,
            "assigned_read_block_rows": 0,
            "unknown_read_block_rows": 0,
            "skipped_regions": [],
        },
        "blocks": [],
    }

    with open_alignment(config) as bam, open_optional_fasta(config.reference_fasta) as fasta:
        references = set(bam.references)
        with gzip.open(marker_path, "wt") as marker_fh, gzip.open(read_tag_path, "wt") as read_fh:
            marker_fh.write("sample_id\tcontig\tpos0\tpos1\tallele1\tallele2\tallele1_count\tallele2_count\tdepth\tthird_allele_fraction\n")
            read_fh.write(
                "sample_id\tread_name\tcontig\tps\tcluster\tmarker_count\tinformative_count\t"
                "score_cluster1\tscore_cluster2\tscore_delta\tconfidence\tstatus\n"
            )
            for requested in regions:
                region = resolve_region(requested, references)
                if region is None:
                    qc["summary"]["skipped_regions"].append({"requested_contig": requested.contig, "reason": "missing_in_alignment_header"})
                    continue
                qc["summary"]["regions_processed"] += 1
                markers = detect_markers(bam, fasta, config, region)
                qc["summary"]["markers_total"] += len(markers)
                for marker in markers:
                    marker_fh.write(
                        f"{config.sample_id}\t{marker.contig}\t{marker.pos0}\t{marker.pos0 + 1}\t"
                        f"{marker.allele1}\t{marker.allele2}\t{marker.allele1_count}\t{marker.allele2_count}\t"
                        f"{marker.depth}\t{marker.third_allele_fraction:.6f}\n"
                    )

                for block_markers in partition_markers(markers, config):
                    block_qc, calls = phase_marker_block(bam, config, block_markers)
                    qc["summary"]["blocks_total"] += 1
                    qc["summary"]["read_block_rows_total"] += len(calls)
                    assigned = sum(1 for call in calls if call.cluster in (1, 2))
                    unknown = len(calls) - assigned
                    qc["summary"]["assigned_read_block_rows"] += assigned
                    qc["summary"]["unknown_read_block_rows"] += unknown
                    if assigned:
                        qc["summary"]["blocks_with_assignments"] += 1
                    qc["blocks"].append(block_qc)
                    for call in calls:
                        read_fh.write(
                            f"{config.sample_id}\t{call.read_name}\t{call.contig}\t{call.ps}\t{call.cluster}\t"
                            f"{call.marker_count}\t{call.informative_count}\t{call.score_cluster1}\t{call.score_cluster2}\t"
                            f"{call.score_delta}\t{call.confidence:.6f}\t{call.status}\n"
                        )

    block_spans = [int(block["span_bp"]) for block in qc["blocks"] if int(block.get("assigned_reads", 0)) > 0]
    qc["summary"]["assigned_block_span_n50"] = n50(block_spans)
    total_rows = int(qc["summary"]["read_block_rows_total"])
    qc["summary"]["assigned_fraction"] = (
        float(qc["summary"]["assigned_read_block_rows"]) / total_rows if total_rows else 0.0
    )
    qc_path.write_text(json.dumps(qc, indent=2, ensure_ascii=False) + "\n")
    return qc


def validate_config(config: BamOnlyPhaseConfig) -> None:
    if config.unknown_read_policy not in {"include_both", "drop_unknown"}:
        raise ValueError("unknown_read_policy must be include_both or drop_unknown")
    if config.min_depth < 1:
        raise ValueError("min_depth must be >= 1")
    if config.min_marker_reads < 1:
        raise ValueError("min_marker_reads must be >= 1")
    if config.min_read_markers < 1:
        raise ValueError("min_read_markers must be >= 1")
    if config.min_score_delta < 0:
        raise ValueError("min_score_delta must be >= 0")
    if not 0.0 <= config.min_minor_allele_fraction <= 0.5:
        raise ValueError("min_minor_allele_fraction must be in [0, 0.5]")
    if not 0.5 <= config.max_major_allele_fraction <= 1.0:
        raise ValueError("max_major_allele_fraction must be in [0.5, 1]")
    if not 0.0 <= config.max_third_allele_fraction <= 1.0:
        raise ValueError("max_third_allele_fraction must be in [0, 1]")
    if config.max_markers_per_block < 2:
        raise ValueError("max_markers_per_block must be >= 2")
    if config.max_marker_gap < 1:
        raise ValueError("max_marker_gap must be >= 1")
    if not config.dry_run:
        if not config.alignment.exists():
            raise FileNotFoundError(config.alignment)
        if config.reference_fasta and not config.reference_fasta.exists():
            raise FileNotFoundError(config.reference_fasta)


def build_contract(config: BamOnlyPhaseConfig) -> dict[str, Any]:
    return {
        "schema": "mfm_longread_bamonly_phase_unoriented_v1",
        "status": "dry_run" if config.dry_run else "planned",
        "sample_id": config.sample_id,
        "alignment": str(config.alignment),
        "output_dir": str(config.output_dir),
        "reference_fasta": str(config.reference_fasta or ""),
        "contigs": config.contigs,
        "regions": config.regions,
        "platform": config.platform,
        "algorithm": {
            "marker_source": "BAM/CRAM count_coverage marker discovery only; no external VCF is used for assignment",
            "marker_type": "biallelic SNV-like A/C/G/T positions only",
            "read_observation_basis": "CIGAR aligned M/=/X query-reference pairs; clipping, insertion, deletion, refskip are not used as marker alleles",
            "cluster_orientation": "arbitrary_unoriented_to_HAP1_HAP2",
            "downstream_policy": "do not treat cluster=1/2 as HP=1/2 until a cluster-to-HAP orientation step passes QC",
        },
        "config": {key: str(value) if isinstance(value, Path) else value for key, value in asdict(config).items()},
    }


def open_alignment(config: BamOnlyPhaseConfig):
    import pysam

    suffix = config.alignment.suffix.lower()
    mode = "rc" if suffix == ".cram" else "rb"
    kwargs: dict[str, Any] = {}
    if config.reference_fasta is not None:
        kwargs["reference_filename"] = str(config.reference_fasta)
    return pysam.AlignmentFile(str(config.alignment), mode, **kwargs)


def open_optional_fasta(path: Path | None):
    if path is None:
        return NullContext(None)
    import pysam

    return pysam.FastaFile(str(path))


class NullContext:
    def __init__(self, value: Any) -> None:
        self.value = value

    def __enter__(self) -> Any:
        return self.value

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> bool:
        return False


def parse_regions(value: str) -> list[Region]:
    regions: list[Region] = []
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            regions.append(Region(normalize_contig(token), None, None))
            continue
        contig_part, interval = token.split(":", 1)
        if "-" not in interval:
            raise ValueError(f"Invalid region interval: {token}")
        left, right = interval.split("-", 1)
        start = max(0, int(left.replace(",", "")) - 1)
        end = int(right.replace(",", ""))
        if end <= start:
            raise ValueError(f"Region end must be greater than start: {token}")
        regions.append(Region(normalize_contig(contig_part), start, end))
    if not regions:
        raise ValueError("regions is empty")
    return regions


def normalize_contig(value: str) -> str:
    value = value.strip()
    if not value:
        raise ValueError("empty contig")
    if value.startswith("chr"):
        return value
    if value in {"X", "Y", "M", "MT"} or value.isdigit():
        return "chrM" if value == "MT" else f"chr{value}"
    return value


def resolve_region(region: Region, references: set[str]) -> Region | None:
    if region.contig in references:
        return region
    no_chr = region.contig.removeprefix("chr")
    if no_chr in references:
        return Region(no_chr, region.start, region.end)
    with_chr = normalize_contig(region.contig)
    if with_chr in references:
        return Region(with_chr, region.start, region.end)
    if region.contig == "chrM" and "MT" in references:
        return Region("MT", region.start, region.end)
    return None


def detect_markers(bam: Any, fasta: Any, config: BamOnlyPhaseConfig, region: Region) -> list[Marker]:
    """Detect candidate SNV-like markers with C-backed base coverage counting.

    The first implementation used ``AlignmentFile.pileup`` and inspected every
    pileup-read object in Python.  That is correct but too slow for real ONT
    CRAM.  ``count_coverage`` computes A/C/G/T coverage in htslib and calls the
    Python read callback once per read, not once per base-column-read pair.  We
    then apply biallelic marker filters vectorially per chunk.
    """
    if config.marker_detection_chunk_bp < 1:
        raise ValueError("marker_detection_chunk_bp must be >= 1")
    import numpy as np

    contig_length = bam.get_reference_length(region.contig)
    start0 = 0 if region.start is None else region.start
    end0 = contig_length if region.end is None else min(region.end, contig_length)
    if end0 <= start0:
        return []

    markers: list[Marker] = []
    base_names = np.array(["A", "C", "G", "T"], dtype=object)

    def read_callback(read: Any) -> bool:
        return is_usable_alignment(read, config)

    chunk_start = start0
    while chunk_start < end0:
        chunk_end = min(chunk_start + config.marker_detection_chunk_bp, end0)
        cov_arrays = bam.count_coverage(
            region.contig,
            start=chunk_start,
            stop=chunk_end,
            quality_threshold=config.min_base_quality,
            read_callback=read_callback,
        )
        counts = np.vstack([np.asarray(arr, dtype=np.int32) for arr in cov_arrays])
        depth = counts.sum(axis=0)
        order = np.argsort(counts, axis=0)
        sorted_counts = np.take_along_axis(counts, order, axis=0)
        major_count = sorted_counts[-1]
        minor_count = sorted_counts[-2]
        third_count = sorted_counts[-3]
        with np.errstate(divide="ignore", invalid="ignore"):
            minor_fraction = np.where(depth > 0, minor_count / depth, 0.0)
            major_fraction = np.where(depth > 0, major_count / depth, 0.0)
            third_fraction = np.where(depth > 0, third_count / depth, 0.0)
        mask = (
            (depth >= config.min_depth)
            & (minor_count >= config.min_marker_reads)
            & (minor_fraction >= config.min_minor_allele_fraction)
            & (major_fraction <= config.max_major_allele_fraction)
            & (third_fraction <= config.max_third_allele_fraction)
        )
        candidate_offsets = np.nonzero(mask)[0]
        for offset in candidate_offsets:
            major_idx = int(order[-1, offset])
            minor_idx = int(order[-2, offset])
            marker = Marker(
                contig=region.contig,
                pos0=chunk_start + int(offset),
                allele1=str(base_names[major_idx]),
                allele2=str(base_names[minor_idx]),
                allele1_count=int(major_count[offset]),
                allele2_count=int(minor_count[offset]),
                depth=int(depth[offset]),
                third_allele_fraction=float(third_fraction[offset]),
            )
            markers.append(marker)
            if config.max_markers_per_contig and len(markers) >= config.max_markers_per_contig:
                return markers
        chunk_start = chunk_end
    return markers


def maybe_build_marker(contig: str, pos0: int, counts: Counter[str], config: BamOnlyPhaseConfig) -> Marker | None:
    depth = sum(counts.values())
    if depth < config.min_depth:
        return None
    if len(counts) < 2:
        return None
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    (allele1, count1), (allele2, count2) = ordered[0], ordered[1]
    if count2 < config.min_marker_reads:
        return None
    minor_fraction = count2 / depth
    major_fraction = count1 / depth
    if minor_fraction < config.min_minor_allele_fraction:
        return None
    if major_fraction > config.max_major_allele_fraction:
        return None
    third_fraction = ordered[2][1] / depth if len(ordered) > 2 else 0.0
    if third_fraction > config.max_third_allele_fraction:
        return None
    return Marker(
        contig=contig,
        pos0=pos0,
        allele1=allele1,
        allele2=allele2,
        allele1_count=count1,
        allele2_count=count2,
        depth=depth,
        third_allele_fraction=third_fraction,
    )


def is_usable_alignment(read: Any, config: BamOnlyPhaseConfig) -> bool:
    if read.is_unmapped or read.is_secondary:
        return False
    if read.is_supplementary and not config.include_supplementary:
        return False
    if config.exclude_duplicates and read.is_duplicate:
        return False
    if read.mapping_quality < config.min_mapq:
        return False
    return True


def partition_markers(markers: list[Marker], config: BamOnlyPhaseConfig) -> Iterable[list[Marker]]:
    current: list[Marker] = []
    previous: Marker | None = None
    for marker in sorted(markers, key=lambda item: item.pos0):
        should_flush = False
        if previous is not None:
            should_flush = (marker.pos0 - previous.pos0 > config.max_marker_gap) or (len(current) >= config.max_markers_per_block)
        if should_flush and current:
            yield current
            current = []
        current.append(marker)
        previous = marker
    if current:
        yield current


def phase_marker_block(bam: Any, config: BamOnlyPhaseConfig, markers: list[Marker]) -> tuple[dict[str, Any], list[ClusterCall]]:
    contig = markers[0].contig
    ps = markers[0].pos0 + 1
    start = markers[0].pos0
    end = markers[-1].pos0 + 1
    allele_index = {marker.pos0: {marker.allele1: 0, marker.allele2: 1} for marker in markers}
    read_vectors = collect_read_vectors(bam, config, contig, start, end, allele_index)
    calls = cluster_read_vectors(config, contig, ps, read_vectors)
    assigned_1 = sum(1 for call in calls if call.cluster == 1)
    assigned_2 = sum(1 for call in calls if call.cluster == 2)
    unknown = sum(1 for call in calls if call.cluster == 0)
    block_qc = {
        "contig": contig,
        "ps": ps,
        "start0": start,
        "end0": end,
        "span_bp": end - start,
        "marker_count": len(markers),
        "reads_with_min_markers": len(read_vectors),
        "assigned_reads": assigned_1 + assigned_2,
        "cluster1_reads": assigned_1,
        "cluster2_reads": assigned_2,
        "unknown_reads": unknown,
        "assigned_fraction": (assigned_1 + assigned_2) / len(calls) if calls else 0.0,
    }
    return block_qc, calls


def collect_read_vectors(
    bam: Any,
    config: BamOnlyPhaseConfig,
    contig: str,
    start: int,
    end: int,
    allele_index: dict[int, dict[str, int]],
) -> dict[str, dict[int, int]]:
    read_vectors: dict[str, dict[int, int]] = defaultdict(dict)
    conflict_positions: dict[str, set[int]] = defaultdict(set)
    marker_positions = set(allele_index)
    for read in bam.fetch(contig, start, end):
        if not is_usable_alignment(read, config) or read.query_sequence is None:
            continue
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            if rpos not in marker_positions or qpos is None:
                continue
            base = read.query_sequence[qpos].upper()
            index = allele_index[rpos].get(base)
            if index is None:
                continue
            previous = read_vectors[read.query_name].get(rpos)
            if previous is not None and previous != index:
                conflict_positions[read.query_name].add(rpos)
                continue
            read_vectors[read.query_name][rpos] = index
    filtered: dict[str, dict[int, int]] = {}
    for read_name, vector in read_vectors.items():
        for conflict_pos in conflict_positions.get(read_name, set()):
            vector.pop(conflict_pos, None)
        if len(vector) >= config.min_read_markers:
            filtered[read_name] = dict(sorted(vector.items()))
    return filtered


def cluster_read_vectors(
    config: BamOnlyPhaseConfig,
    contig: str,
    ps: int,
    read_vectors: dict[str, dict[int, int]],
) -> list[ClusterCall]:
    if not read_vectors:
        return []
    if len(read_vectors) == 1:
        read_name, vector = next(iter(read_vectors.items()))
        return [
            ClusterCall(
                read_name=read_name,
                contig=contig,
                ps=ps,
                cluster=0,
                marker_count=len(vector),
                informative_count=0,
                score_cluster1=0,
                score_cluster2=0,
                score_delta=0,
                confidence=0.0,
                status="unknown_single_read_block",
            )
        ]

    profile1, profile2 = initialize_profiles(read_vectors)
    assignments: dict[str, tuple[int, int, int, int, int, str]] = {}
    for _ in range(config.max_iterations):
        assignments = assign_reads(config, read_vectors, profile1, profile2)
        new_profile1, new_profile2 = update_profiles(read_vectors, assignments, profile1, profile2)
        if new_profile1 == profile1 and new_profile2 == profile2:
            break
        profile1, profile2 = new_profile1, new_profile2

    calls: list[ClusterCall] = []
    for read_name in sorted(read_vectors):
        vector = read_vectors[read_name]
        cluster, score1, score2, delta, informative, status = assignments.get(
            read_name, (0, 0, 0, 0, 0, "unknown_unassigned")
        )
        calls.append(
            ClusterCall(
                read_name=read_name,
                contig=contig,
                ps=ps,
                cluster=cluster,
                marker_count=len(vector),
                informative_count=informative,
                score_cluster1=score1,
                score_cluster2=score2,
                score_delta=delta,
                confidence=(delta / informative) if informative else 0.0,
                status=status,
            )
        )
    return calls


def initialize_profiles(read_vectors: dict[str, dict[int, int]]) -> tuple[dict[int, int], dict[int, int]]:
    seed_name = max(read_vectors, key=lambda name: (len(read_vectors[name]), name))
    seed_vector = read_vectors[seed_name]
    profile1 = dict(seed_vector)
    profile2 = {pos: 1 - allele for pos, allele in seed_vector.items()}

    best_opposite: tuple[int, str] | None = None
    for read_name, vector in read_vectors.items():
        if read_name == seed_name:
            continue
        shared = set(seed_vector).intersection(vector)
        if not shared:
            continue
        discordant = sum(1 for pos in shared if seed_vector[pos] != vector[pos])
        concordant = len(shared) - discordant
        score = discordant - concordant
        if best_opposite is None or score > best_opposite[0]:
            best_opposite = (score, read_name)
    if best_opposite and best_opposite[0] > 0:
        for pos, allele in read_vectors[best_opposite[1]].items():
            profile2[pos] = allele
            profile1.setdefault(pos, 1 - allele)
    return profile1, profile2


def assign_reads(
    config: BamOnlyPhaseConfig,
    read_vectors: dict[str, dict[int, int]],
    profile1: dict[int, int],
    profile2: dict[int, int],
) -> dict[str, tuple[int, int, int, int, int, str]]:
    assignments: dict[str, tuple[int, int, int, int, int, str]] = {}
    for read_name, vector in read_vectors.items():
        score1, info1 = score_against_profile(vector, profile1)
        score2, info2 = score_against_profile(vector, profile2)
        informative = max(info1, info2)
        delta = abs(score1 - score2)
        if informative < config.min_read_markers:
            assignments[read_name] = (0, score1, score2, delta, informative, "unknown_insufficient_informative_markers")
            continue
        if delta < config.min_score_delta or score1 == score2:
            assignments[read_name] = (0, score1, score2, delta, informative, "unknown_low_score_delta")
            continue
        cluster = 1 if score1 > score2 else 2
        assignments[read_name] = (cluster, score1, score2, delta, informative, "assigned")
    return assignments


def score_against_profile(vector: dict[int, int], profile: dict[int, int]) -> tuple[int, int]:
    score = 0
    informative = 0
    for pos, allele in vector.items():
        expected = profile.get(pos)
        if expected is None:
            continue
        informative += 1
        score += 1 if allele == expected else -1
    return score, informative


def update_profiles(
    read_vectors: dict[str, dict[int, int]],
    assignments: dict[str, tuple[int, int, int, int, int, str]],
    previous1: dict[int, int],
    previous2: dict[int, int],
) -> tuple[dict[int, int], dict[int, int]]:
    counts: dict[int, dict[int, Counter[int]]] = {1: defaultdict(Counter), 2: defaultdict(Counter)}
    for read_name, vector in read_vectors.items():
        cluster = assignments.get(read_name, (0, 0, 0, 0, 0, ""))[0]
        if cluster not in (1, 2):
            continue
        for pos, allele in vector.items():
            counts[cluster][pos][allele] += 1
    all_positions = set(previous1) | set(previous2) | set(counts[1]) | set(counts[2])
    profile1: dict[int, int] = {}
    profile2: dict[int, int] = {}
    for pos in all_positions:
        allele1 = majority_allele(counts[1].get(pos))
        allele2 = majority_allele(counts[2].get(pos))
        if allele1 is None and allele2 is None:
            allele1 = previous1.get(pos)
            allele2 = previous2.get(pos)
        if allele1 is None and allele2 is not None:
            allele1 = 1 - allele2
        if allele2 is None and allele1 is not None:
            allele2 = 1 - allele1
        if allele1 is None or allele2 is None:
            continue
        if allele1 == allele2:
            old1 = previous1.get(pos)
            old2 = previous2.get(pos)
            if old1 is not None and old2 is not None and old1 != old2:
                allele1, allele2 = old1, old2
            else:
                allele2 = 1 - allele1
        profile1[pos] = allele1
        profile2[pos] = allele2
    return profile1, profile2


def majority_allele(counter: Counter[int] | None) -> int | None:
    if not counter:
        return None
    return sorted(counter.items(), key=lambda item: (-item[1], item[0]))[0][0]


def n50(lengths: Iterable[int]) -> int:
    ordered = sorted((length for length in lengths if length > 0), reverse=True)
    if not ordered:
        return 0
    half = sum(ordered) / 2
    cumulative = 0
    for length in ordered:
        cumulative += length
        if cumulative >= half:
            return length
    return ordered[-1]


# Lightweight self-test entry point used by validation smoke commands.  It is
# intentionally kept dependency-free so py_compile/unit checks do not require a
# real BAM fixture.
def _synthetic_cluster_smoke() -> dict[str, Any]:
    config = BamOnlyPhaseConfig(
        sample_id="synthetic",
        alignment=Path("/tmp/synthetic.bam"),
        output_dir=Path("/tmp"),
        min_read_markers=2,
        min_score_delta=2,
    )
    read_vectors = {
        "r1": {10: 0, 20: 0, 30: 1, 40: 0},
        "r2": {10: 0, 20: 0, 30: 1, 40: 0},
        "r3": {10: 1, 20: 1, 30: 0, 40: 1},
        "r4": {10: 1, 20: 1, 30: 0, 40: 1},
        "r5": {10: 0, 20: 1},
    }
    calls = cluster_read_vectors(config, "chrS", 11, read_vectors)
    assigned = {call.read_name: call.cluster for call in calls}
    return {
        "assigned": assigned,
        "cluster_counts": dict(Counter(assigned.values())),
        "pass": assigned["r1"] == assigned["r2"] and assigned["r3"] == assigned["r4"] and assigned["r1"] != assigned["r3"],
    }
