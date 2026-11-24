from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Any

import pysam

from .config import (
    MIN_ILLUMINA_PAIRED_PRIMARY_FRAC,
    PersonalRealignConfig,
    PreparedPairedInput,
)
from .policy import is_hap_specific_mode, is_illumina_platform

def illumina_prepared_path(work_dir: Path, sample: str, route_label: str) -> Path:
    return work_dir / f"{sample}.illumina.{route_label}.paired_source.bam"
def primary_hp(read: pysam.AlignedSegment) -> int:
    if not read.has_tag("HP"):
        return 0
    try:
        value = int(read.get_tag("HP"))
    except (TypeError, ValueError):
        return -1
    return value if value in {0, 1, 2} else -1


def classify_illumina_fragment(reads: list[pysam.AlignedSegment]) -> tuple[str, int]:
    """Return a mutually-exclusive QNAME classification for strict routing."""
    mapped_reads = [read for read in reads if not read.is_unmapped]
    if not mapped_reads:
        return "unknown", 0
    hps = {primary_hp(read) for read in mapped_reads}
    if -1 in hps:
        return "invalid_hp", 0
    tagged = hps & {1, 2}
    if tagged == {1}:
        return "hap1", 1
    if tagged == {2}:
        return "hap2", 2
    if tagged == {1, 2}:
        return "conflict", 0
    return "unknown", 0


def iter_qname_groups(source: pysam.AlignmentFile):
    current_name: str | None = None
    group: list[pysam.AlignedSegment] = []
    for read in source.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue
        if current_name is None:
            current_name = read.query_name
        if read.query_name != current_name:
            if group:
                yield group
            current_name = read.query_name
            group = []
        group.append(read)
    if group:
        yield group


def carry_source_duplicate(
    config: PersonalRealignConfig, read: pysam.AlignedSegment, stats: dict[str, int]
) -> None:
    if not read.is_duplicate:
        return
    stats["source_duplicate_primary"] = stats.get("source_duplicate_primary", 0) + 1
    if config.illumina_duplicate_mode == "preserve_source":
        read.set_tag("OD", 1, value_type="i", replace=True)


def add_alignment_input_metrics(path: Path, stats: dict[str, Any], threads: int = 1, *,
                                workers: int | None = None, samtools_bin: str = "samtools",
                                chunk_bytes: int = 8 * 1024 ** 2, concurrency: int = 1) -> None:
    """소형 입력은 직렬, 큰 입력은 QNAME 경계를 보존하는 병렬 스트림으로 집계함."""
    import json
    import resource
    import time
    from .input_metrics import metric_plan, parallel_input_metrics
    plan = metric_plan(path, threads, workers=workers, chunk_bytes=chunk_bytes, concurrency=concurrency)
    print("[INPUT-METRICS] " + json.dumps({"phase": "plan", "path": str(path), **plan}), flush=True)
    started, cpu_started = time.monotonic(), time.process_time()
    if plan["mode"] == "parallel_sam_stream":
        metrics, execution = parallel_input_metrics(path, samtools_bin, plan)
        stats.update(metrics)
    else:
        _add_alignment_input_metrics_serial(path, stats, threads=plan["samtools_aux_threads"])
        execution = {**plan, "workers_used": 1, "batches": 1,
                     "wall_seconds": round(time.monotonic() - started, 3),
                     "worker_cpu_seconds": round(time.process_time() - cpu_started, 3),
                     "worker_peak_rss_sum_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024}
    stats["input_metrics_execution"] = execution
    print("[INPUT-METRICS] " + json.dumps({"phase": "complete", **execution}), flush=True)


def _add_alignment_input_metrics_serial(path: Path, stats: dict[str, int], threads: int = 1) -> None:
    """Count exactly the complete R1/R2 records emitted by samtools fastq."""
    metrics = {
        "paired_fragments": 0,
        "alignment_input_primary": 0,
        "alignment_input_routed_primary": 0,
        "alignment_input_rescue_primary": 0,
        "alignment_input_hp1": 0,
        "alignment_input_hp2": 0,
        "alignment_input_untagged": 0,
        "alignment_input_duplicate": 0,
        "singleton_discard_records": 0,
        "read_other_discard_records": 0,
    }
    with pysam.AlignmentFile(str(path), "rb", threads=max(0, threads)) as source:
        for reads in iter_qname_groups(source):
            read1 = [read for read in reads if read.is_read1 and not read.is_read2]
            read2 = [read for read in reads if read.is_read2 and not read.is_read1]
            other = [read for read in reads if read not in read1 and read not in read2]
            if read1 and read2:
                selected = (read1[0], read2[0])
                metrics["paired_fragments"] += 1
                metrics["alignment_input_primary"] += 2
                for read in selected:
                    if read.has_tag("XU") and int(read.get_tag("XU")) == 1:
                        metrics["alignment_input_rescue_primary"] += 1
                    else:
                        metrics["alignment_input_routed_primary"] += 1
                        hp = primary_hp(read)
                        if hp in {1, 2}:
                            metrics[f"alignment_input_hp{hp}"] += 1
                        else:
                            metrics["alignment_input_untagged"] += 1
                    metrics["alignment_input_duplicate"] += int(read.has_tag("OD") and int(read.get_tag("OD")) == 1)
                metrics["singleton_discard_records"] += len(read1[1:]) + len(read2[1:])
            else:
                metrics["singleton_discard_records"] += len(read1) + len(read2)
            metrics["read_other_discard_records"] += len(other)
    stats.update(metrics)


def prepare_illumina_hap_specific_inputs(
    config: PersonalRealignConfig,
    sample: str,
    alignment: Path,
    contigs: list[str],
    work_dir: Path,
) -> tuple[dict[int, PreparedPairedInput], dict[str, Any]]:
    """Route Illumina QNAME fragments to exactly one HAP-specific source BAM."""
    outputs = {hap: illumina_prepared_path(work_dir, sample, f"hap{hap}_specific") for hap in (1, 2)}
    parts = {hap: path.with_suffix(path.suffix + ".part") for hap, path in outputs.items()}
    name_sorted = work_dir / f"{sample}.illumina.hap_specific.name_sorted.bam"
    sort_prefix = work_dir / f"{sample}.illumina.hap_specific.name_sort"
    for path in [*parts.values(), name_sorted]:
        if path.exists():
            path.unlink()
    command = [
        config.samtools_bin, "sort", "-N", "-@", str(max(1, config.sort_threads)),
        "-m", config.sort_mem, "-T", str(sort_prefix), "-o", str(name_sorted),
    ]
    if alignment.suffix.lower() == ".cram" and config.reference_fasta is not None:
        command.extend(["--reference", str(config.reference_fasta)])
    command.append(str(alignment))
    subprocess.run(command, check=True)

    stats_by_hap = {
        hap: {
            "routed_primary": 0, "rescue_unmapped_primary": 0,
            "routed_hp1": 0, "routed_hp2": 0, "routed_untagged": 0,
            "retained_fragments": 0,
            "source_duplicate_primary": 0,
        }
        for hap in (1, 2)
    }
    qc: dict[str, Any] = {
        "mode": "hap_specific",
        "granularity": "QNAME_fragment",
        "retention_rule": "HP1-only -> HAP1; HP2-only -> HAP2; HP0/missing, HP1/HP2 conflict, invalid HP -> drop",
        "source_fragments": 0,
        "hap1_fragments": 0,
        "hap2_fragments": 0,
        "unknown_fragments": 0,
        "conflict_fragments": 0,
        "invalid_hp_fragments": 0,
        "one_mate_hp_evidence_fragments": 0,
        "unmapped_without_hp_evidence_fragments": 0,
        "qname_intersection_hap1_hap2": 0,
        "qname_intersection_assertion": "PASS: each QNAME dispatches through one mutually exclusive branch",
    }
    contig_set = set(contigs)
    try:
        with pysam.AlignmentFile(
            str(name_sorted), "rb", threads=htslib_decode_threads(config)
        ) as source:
            with pysam.AlignmentFile(str(parts[1]), "wb", header=source.header) as dest1, pysam.AlignmentFile(
                str(parts[2]), "wb", header=source.header
            ) as dest2:
                destinations = {1: dest1, 2: dest2}
                for reads in iter_qname_groups(source):
                    qc["source_fragments"] += 1
                    label, hap = classify_illumina_fragment(reads)
                    if label not in {"hap1", "hap2"}:
                        qc[f"{label}_fragments"] += 1
                        if label == "unknown" and all(read.is_unmapped for read in reads):
                            qc["unmapped_without_hp_evidence_fragments"] += 1
                        continue
                    if sum(primary_hp(read) in {1, 2} for read in reads) == 1:
                        qc["one_mate_hp_evidence_fragments"] += 1
                    eligible_mapped = [
                        read for read in reads if not read.is_unmapped and read.reference_name in contig_set
                    ]
                    if not eligible_mapped:
                        continue
                    eligible = [*eligible_mapped, *(read for read in reads if read.is_unmapped)]
                    qc[f"hap{hap}_fragments"] += 1
                    stats = stats_by_hap[hap]
                    stats["retained_fragments"] += 1
                    for read in eligible:
                        if read.is_unmapped:
                            read.set_tag("XU", 1, value_type="i", replace=True)
                            stats["rescue_unmapped_primary"] += 1
                        else:
                            stats["routed_primary"] += 1
                            hp = primary_hp(read)
                            if hp == 1:
                                stats["routed_hp1"] += 1
                            elif hp == 2:
                                stats["routed_hp2"] += 1
                            else:
                                stats["routed_untagged"] += 1
                        carry_source_duplicate(config, read, stats)
                        destinations[hap].write(read)
    except Exception:
        for path in parts.values():
            if path.exists():
                path.unlink()
        raise
    finally:
        if name_sorted.exists():
            name_sorted.unlink()
    for hap in (1, 2):
        parts[hap].replace(outputs[hap])
        add_alignment_input_metrics(outputs[hap], stats_by_hap[hap], threads=htslib_decode_threads(config),
                                    samtools_bin=config.samtools_bin, concurrency=config.jobs)
    qname_intersection = prepared_qname_intersection_count(outputs[1], outputs[2])
    qc["qname_intersection_hap1_hap2"] = qname_intersection
    if qname_intersection:
        for path in outputs.values():
            if path.exists():
                path.unlink()
        raise RuntimeError(f"strict Illumina routing produced {qname_intersection} overlapping HAP QNAME(s)")
    prepared = {
        hap: PreparedPairedInput(
            path=outputs[hap], route_label=f"hap{hap}_specific",
            stats=stats_by_hap[hap], name_collated=True,
        )
        for hap in (1, 2)
    }
    qc["per_haplotype"] = {f"HAP{hap}": stats_by_hap[hap] for hap in (1, 2)}
    return prepared, qc


def prepared_qname_intersection_count(hap1_path: Path, hap2_path: Path) -> int:
    """Count QNAME overlap between two name-collated strict prepared BAMs."""
    def names(path: Path):
        previous: str | None = None
        with pysam.AlignmentFile(str(path), "rb") as source:
            for read in source.fetch(until_eof=True):
                if read.is_secondary or read.is_supplementary or read.query_name == previous:
                    continue
                previous = read.query_name
                yield previous

    left = iter(names(hap1_path))
    right = iter(names(hap2_path))
    value1 = next(left, None)
    value2 = next(right, None)
    overlap = 0
    while value1 is not None and value2 is not None:
        if value1 == value2:
            overlap += 1
            value1 = next(left, None)
            value2 = next(right, None)
        elif value1 < value2:
            value1 = next(left, None)
        else:
            value2 = next(right, None)
    return overlap


def select_read_for_haplotype(config: PersonalRealignConfig, read: pysam.AlignedSegment, hap: int) -> bool:
    if config.read_routing_policy == "all_to_both":
        return True
    if read.has_tag("HP") and int(read.get_tag("HP")) == hap:
        return True
    return config.unknown_read_policy == "include_both" and (not read.has_tag("HP") or int(read.get_tag("HP")) == 0)


def htslib_decode_threads(config: PersonalRealignConfig) -> int:
    # BAM/CRAM block decompression만 병렬화한다. record 선택과 QC 집계 순서는
    # 기존 단일 Python iterator를 유지하므로 read/QC 결과는 변하지 않는다.
    return max(0, min(4, int(config.sort_threads)))


def empty_input_counts() -> dict[str, int]:
    return {"hp1": 0, "hp2": 0, "untagged": 0, "all_primary_mapped": 0, "unmapped": 0}


def count_input_reads(config: PersonalRealignConfig, alignment: Path, contigs: list[str]) -> dict[str, int]:
    if not config.verify_counts:
        return empty_input_counts()
    if is_illumina_platform(config.platform):
        return count_illumina_input_reads(config, alignment, contigs)
    return count_longread_input_reads(config, alignment, contigs)


def count_longread_input_reads(
    config: PersonalRealignConfig, alignment: Path, contigs: list[str]
) -> dict[str, int]:
    """Long-read HP classes from one pass, matching the old samtools expressions."""
    return _count_input_reads_one_pass(config, alignment, contigs, invalid_hp_as_untagged=False)


def count_illumina_input_reads(
    config: PersonalRealignConfig, alignment: Path, contigs: list[str]
) -> dict[str, int]:
    """Illumina routing counts from one pass, preserving its prior unknown-HP rule."""
    return _count_input_reads_one_pass(config, alignment, contigs, invalid_hp_as_untagged=True)


def _count_input_reads_one_pass(
    config: PersonalRealignConfig,
    alignment: Path,
    contigs: list[str],
    *,
    invalid_hp_as_untagged: bool,
    _force_serial: bool = False,
) -> dict[str, int]:
    import os
    from ..input_parallel import cpu_budget
    backend=os.environ.get('MFM_REALIGN_INPUT_BACKEND','auto')
    if backend not in {'auto','serial','indexed'}:raise ValueError('MFM_REALIGN_INPUT_BACKEND는 auto|serial|indexed여야함')
    live_count_cpu=cpu_budget()
    if os.environ.get('MFM_CPU_POOL_LEASE_FILE'):
        import json
        live_count_cpu=int(json.loads(Path(os.environ['MFM_CPU_POOL_LEASE_FILE']).read_text())['cpus'])
    if not _force_serial and (backend=='indexed' or (backend=='auto' and getattr(config,'jobs',1)==1 and alignment.stat().st_size>=32*1024**2 and live_count_cpu>1)):
        if getattr(config,'jobs',1)>1 and os.environ.get('MFM_CPU_POOL_LEASE_FILE'):raise ValueError('공유lease의indexed입력집계는단일sample호출에서만지원함')
        with pysam.AlignmentFile(str(alignment),reference_filename=str(config.reference_fasta) if config.reference_fasta else None) as handle:indexed=handle.has_index()
        if indexed:
            from .input_count import count_indexed
            return count_indexed(config,alignment,contigs,invalid_hp_as_untagged)
        if backend=='indexed':raise ValueError('indexed입력집계에BAM/CRAM index가필요함')
    counts = empty_input_counts()
    contig_set = set(contigs)
    open_kwargs: dict[str, Any] = {"threads": htslib_decode_threads(config)}
    if alignment.suffix.lower() == ".cram" and config.reference_fasta is not None:
        open_kwargs["reference_filename"] = str(config.reference_fasta)
    with pysam.AlignmentFile(str(alignment), "rb", **open_kwargs) as source:
        for read in source.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            if read.is_unmapped:
                counts["unmapped"] += 1
                continue
            if read.reference_name not in contig_set:
                continue
            counts["all_primary_mapped"] += 1
            if not read.has_tag("HP"):
                counts["untagged"] += 1
                continue
            try:
                hp = int(read.get_tag("HP"))
            except (TypeError, ValueError):
                if invalid_hp_as_untagged:
                    raise
                continue
            if hp == 1:
                counts["hp1"] += 1
            elif hp == 2:
                counts["hp2"] += 1
            elif hp == 0 or invalid_hp_as_untagged:
                counts["untagged"] += 1
    return counts
def samtools_count(config: PersonalRealignConfig, alignment: Path, flag_opt: str, flag_value: str, expr_opt: str | None, expr: str | None, regions: list[str]) -> int:
    cmd = [config.samtools_bin, "view", "-c", flag_opt, flag_value]
    if expr_opt and expr:
        cmd.extend([expr_opt, expr])
    if config.reference_fasta is not None:
        cmd.extend(["-T", str(config.reference_fasta)])
    cmd.append(str(alignment))
    cmd.extend(regions)
    result = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return int(result.stdout.strip() or 0)


def scan_realigned_output(config, output, target_fa, hap, *, workers=1, execution=None, resource_limits=None):
    from .output_scan import scan_realigned_output as scan
    return scan(config,output,target_fa,hap,workers=workers,execution=execution,resource_limits=resource_limits)


from .output_scan import format_flagstat_text


def write_realign_counts(path: Path, scan: dict[str, Any], hap: int) -> None:
    rows = [
        ("primary", scan["primary"]),
        ("secondary", scan["secondary"]),
        ("supplementary", scan["supplementary"]),
        ("paired_primary", scan["paired_primary"]),
        ("proper_pair_primary", scan["proper_pair_primary"]),
        ("paired_primary_frac", scan["paired_primary_frac"]),
        ("duplicate_primary", scan["duplicate_primary"]),
        ("od_tag_primary", scan["od_tag_primary"]),
        ("routed_primary", scan["routed_primary"]),
        ("routed_mapped", scan["routed_mapped"]),
        ("routed_hp1", scan["routed_hp1"]),
        ("routed_hp2", scan["routed_hp2"]),
        ("routed_untagged", scan["routed_untagged"]),
        ("rescue_primary", scan["rescue_primary"]),
        ("rescue_mapped", scan["rescue_mapped"]),
        ("routed_mapped_frac", scan["routed_mapped_frac"]),
        ("rescue_mapped_frac", scan["rescue_mapped_frac"]),
        ("nm_sum", scan["nm_sum"]),
        ("mapped_bp", scan["mapped_bp"]),
        ("nm_per_base_all_mapped", scan["nm_per_base_all_mapped"]),
    ]
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def expected_routed_primary(config: PersonalRealignConfig, hap: int, input_counts: dict[str, int]) -> int:
    if config.read_routing_policy == "all_to_both":
        return input_counts["all_primary_mapped"]
    if config.unknown_read_policy == "include_both":
        return input_counts[f"hp{hap}"] + input_counts["untagged"]
    return input_counts[f"hp{hap}"]


def verify_realign_counts(
    config: PersonalRealignConfig,
    hap: int,
    input_counts: dict[str, int],
    scan: dict[str, Any],
    prepared_input: PreparedPairedInput | None,
) -> None:
    require_illumina_pairs = prepared_input is not None
    if config.verify_counts and prepared_input is not None:
        stats = prepared_input.stats
        require_count_within_loss_tolerance(
            f"HAP{hap} paired-input routed primary",
            scan["routed_primary"],
            stats["alignment_input_routed_primary"],
            config.count_tolerance_abs,
        )
        require_count_within_loss_tolerance(
            f"HAP{hap} paired-input XU rescue primary",
            scan["rescue_primary"],
            stats["alignment_input_rescue_primary"],
            config.count_tolerance_abs,
        )
        for hp in (1, 2):
            require_count_within_loss_tolerance(
                f"HAP{hap} paired-input HP{hp}", scan[f"routed_hp{hp}"],
                stats[f"alignment_input_hp{hp}"], config.count_tolerance_abs,
            )
        require_count_within_loss_tolerance(
            f"HAP{hap} paired-input untagged", scan["routed_untagged"],
            stats["alignment_input_untagged"], config.count_tolerance_abs,
        )
        expected_duplicates = stats["alignment_input_duplicate"] if config.illumina_duplicate_mode == "preserve_source" else 0
        require_count_within_loss_tolerance(
            f"HAP{hap} source duplicate flag", scan["duplicate_primary"], expected_duplicates, 0
        )
        if scan["od_tag_primary"]:
            raise RuntimeError(f"HAP{hap} output retains {scan['od_tag_primary']} internal OD tag(s)")
        if is_hap_specific_mode(config) and scan[f"routed_hp{3 - hap}"] != 0:
            raise RuntimeError(f"HAP{hap} strict output contains HP:i:{3 - hap} routed reads")
    elif not config.verify_counts:
        if scan[f"routed_hp{hap}"] <= 0:
            raise RuntimeError(f"HAP{hap} output has no preserved HP:i:{hap} reads")
    else:
        expected_routed = expected_routed_primary(config, hap, input_counts)
        require_count_within_loss_tolerance(
            f"HAP{hap} routed primary", scan["routed_primary"], expected_routed, config.count_tolerance_abs
        )
        require_count_within_loss_tolerance(
            f"HAP{hap} XU rescue primary",
            scan["rescue_primary"],
            input_counts["unmapped"] if config.include_unmapped else 0,
            config.count_tolerance_abs,
        )
        if config.read_routing_policy == "hp_routed":
            require_count_within_loss_tolerance(
                f"HAP{hap} HP tag count", scan[f"routed_hp{hap}"], input_counts[f"hp{hap}"], config.count_tolerance_abs
            )
        if config.read_routing_policy == "hp_routed":
            expected_unknown = input_counts["untagged"] if config.unknown_read_policy == "include_both" else 0
            require_count_within_loss_tolerance(
                f"HAP{hap} untagged count under unknown_read_policy={config.unknown_read_policy}",
                scan["routed_untagged"],
                expected_unknown,
                config.count_tolerance_abs,
            )
        if config.read_routing_policy == "all_to_both":
            for hp in (1, 2):
                require_count_within_loss_tolerance(
                    f"HAP{hap} HP{hp} count", scan[f"routed_hp{hp}"], input_counts[f"hp{hp}"], config.count_tolerance_abs
                )
            require_count_within_loss_tolerance(
                f"HAP{hap} untagged count", scan["routed_untagged"], input_counts["untagged"], config.count_tolerance_abs
            )
    if config.ratio_qc_policy == "strict":
        if scan["routed_mapped_frac"] < config.min_mapped_frac:
            raise RuntimeError(f"HAP{hap} routed mapped fraction {scan['routed_mapped_frac']} < {config.min_mapped_frac}")
        if require_illumina_pairs and scan["paired_primary_frac"] < MIN_ILLUMINA_PAIRED_PRIMARY_FRAC:
            raise RuntimeError(
                f"HAP{hap} paired primary fraction {scan['paired_primary_frac']} < {MIN_ILLUMINA_PAIRED_PRIMARY_FRAC}"
            )


def require_count_within_loss_tolerance(label: str, actual: int, expected: int, tolerance: int) -> None:
    lost = expected - actual
    if lost < 0 or lost > tolerance:
        raise RuntimeError(f"{label} {actual} != input {expected} (lost={lost}, allowed=0..{tolerance})")


def nm_per_base(config: PersonalRealignConfig, path: Path, reference_fasta: Path | None, contig: str, max_reads: int) -> dict[str, Any]:
    kwargs: dict[str, Any] = {"threads": htslib_decode_threads(config)}
    if reference_fasta is not None:
        kwargs["reference_filename"] = str(reference_fasta)
    nm = 0
    bp = 0
    reads = 0
    with pysam.AlignmentFile(str(path), "rb", **kwargs) as bam:
        length = bam.get_reference_length(contig)
        start, end = (20_000_000, 25_000_000) if length and length > 25_000_000 else (0, length)
        for read in bam.fetch(contig, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.has_tag("NM"):
                nm += int(read.get_tag("NM"))
            bp += int(read.query_length or len(read.query_sequence or ""))
            reads += 1
            if reads >= max_reads:
                break
    return {"nm_per_base": round(nm / max(bp, 1), 8), "reads": reads, "bp": bp, "nm": nm}
