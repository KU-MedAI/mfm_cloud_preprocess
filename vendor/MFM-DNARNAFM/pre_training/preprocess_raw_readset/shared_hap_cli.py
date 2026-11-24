"""Shared-read HAP1/HAP2 chromosome preprocessing CLI.

Input TSV columns:
  sample, contig, hap1_fasta, hap1_chain, hap2_fasta, hap2_chain, bam

Output layout:
  <output-root>/<sample>/<contig>.zarr
    reads/
    HAP1/
    HAP2/
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .config import ExportConfig
from .shared_hap.manifest import (
    init_manifest,
    mark_manifest_complete,
    mark_manifest_failed,
    update_manifest_record,
    write_json,
)
from .shared_hap.records import hap_paths, read_record_tsv
from .shared_hap_writer import (
    HAP_LABELS,
    SharedHapChrWriter,
    build_paired_hap_window_artifacts,
    prepare_hap_context,
    validate_shared_hap_zarr,
)

# GRCh38 primary chromosome sizes. Used only for job scheduling; output
# coordinates and windows are still derived from each HAP FASTA/chain.
_GRCH38_CONTIG_SIZE = {
    "chr1": 248_956_422, "chr2": 242_193_529, "chr3": 198_295_559,
    "chr4": 190_214_555, "chr5": 181_538_259, "chr6": 170_805_979,
    "chr7": 159_345_973, "chr8": 145_138_636, "chr9": 138_394_717,
    "chr10": 133_797_422, "chr11": 135_086_622, "chr12": 133_275_309,
    "chr13": 114_364_328, "chr14": 107_043_718, "chr15": 101_991_189,
    "chr16": 90_338_345, "chr17": 83_257_441, "chr18": 80_373_285,
    "chr19": 58_617_616, "chr20": 64_444_167, "chr21": 46_709_983,
    "chr22": 50_818_468,
}


def _job_schedule_key(job: dict) -> tuple[int, str, str]:
    """Sort larger contigs first to reduce the multiprocessing tail.

    With many samples, submitting jobs in TSV order can leave a few large
    chromosomes at the end if the TSV is not already balanced. Longest-first
    scheduling keeps workers occupied more evenly. This does not change any
    output content; it only changes execution order.
    """
    contig = str(job["contig"])
    normalized = contig if contig.startswith("chr") else f"chr{contig}"
    size = _GRCH38_CONTIG_SIZE.get(normalized, 0)
    return (-int(size), str(job["sample"]), contig)



def _existing_pass(path: Path, window_size: int, validation_level: str) -> bool:
    if not path.exists() or not (path / "chr_qc.json").exists():
        return False
    try:
        qc = validate_shared_hap_zarr(path, window_size=window_size, validation_level=validation_level)
    except Exception:
        return False
    return qc.get("status") == "pass"



def _run_contig(job: dict) -> dict:
    started = time.time()
    config = ExportConfig(**job["config"])
    sample = job["sample"]
    contig = job["contig"]
    out_path = Path(job["output_root"]) / sample / f"{contig}.zarr"
    partial = out_path.parent / f".{contig}.zarr.partial.{os.getpid()}"
    if partial.exists():
        shutil.rmtree(partial)
    if out_path.exists():
        if job["overwrite"]:
            shutil.rmtree(out_path)
        elif job["resume"] and _existing_pass(out_path, config.window_size, job["validation_level"]):
            qc = json.loads((out_path / "chr_qc.json").read_text(encoding="utf-8"))
            return {
                "status": "skip_qc_pass",
                "sample": sample,
                "contig": contig,
                "zarr": str(out_path),
                "shared_read_count": qc.get("shared_read_count", 0),
                "haplotypes": qc.get("haplotypes", {}),
                "read_haplotype_counts": qc.get("read_haplotype_counts", {}),
                "read_phase_set": qc.get("read_phase_set", {}),
                "elapsed_seconds": 0.0,
            }
        else:
            raise FileExistsError(f"existing output zarr: {out_path}")

    try:
        summaries: dict[str, dict] = {}
        hap_window_counts: dict[str, int] = {}
        hap_contexts: dict[str, dict] = {}
        for label in HAP_LABELS:
            fasta, chain = hap_paths(job, label)
            context = prepare_hap_context(
                sample=sample,
                label=label,
                contig=contig,
                individual_fasta=fasta,
                chain_path=chain,
                config=config,
            )
            total = int(context["summary"]["coordinate_eligible_window_count_total"])
            if job.get("max_windows") is not None:
                total = min(total, int(job["max_windows"]))
            if total <= 0:
                raise ValueError(f"no eligible windows for {sample}/{label}/{contig}")
            summaries[label] = context["summary"]
            hap_window_counts[label] = total
            hap_contexts[label] = context

        writer = SharedHapChrWriter(
            partial,
            sample=sample,
            contig=contig,
            hap_window_counts=hap_window_counts,
            window_size=config.window_size,
            attrs={
                "source_bam_path": str(job["bam"]),
                "reference_fasta": str(job["reference_fasta"]) if job.get("reference_fasta") else None,
                "source_hap1_fasta": str(job["hap1_fasta"]),
                "source_hap1_chain": str(job["hap1_chain"]),
                "source_hap2_fasta": str(job["hap2_fasta"]),
                "source_hap2_chain": str(job["hap2_chain"]),
                "config": config.to_dict(),
                "window_summary": summaries,
                "builder": "paired_hap_union_bam_fetch_hp_ps",
            },
        )
        chunk_size = int(job["windows_per_chunk"])
        max_total = max(hap_window_counts.values())
        for start in range(0, max_total, chunk_size):
            artifacts_by_label = build_paired_hap_window_artifacts(
                sample=sample,
                contig=contig,
                bam_path=Path(job["bam"]),
                hap_contexts=hap_contexts,
                reference_fasta=Path(job["reference_fasta"]) if job.get("reference_fasta") else None,
                max_windows=min(chunk_size, max_total - start),
                window_start_index=start,
                config=config,
            )
            for label, artifacts in artifacts_by_label.items():
                # `max_windows` may cap one HAP earlier than the other if their
                # eligible window counts differ. The paired builder simply omits
                # labels without selected windows for this chunk.
                if writer.hap_window_cursor[label] < hap_window_counts[label]:
                    writer.append_hap_chunk(label, artifacts)
        qc = writer.finalize(validation_level=job["validation_level"])
        if qc["status"] != "pass":
            raise RuntimeError(f"shared HAP chr QC failed for {sample}/{contig}: {qc['errors']}")
        partial.rename(out_path)
        qc["zarr"] = str(out_path)
        qc["elapsed_seconds"] = round(time.time() - started, 3)
        write_json(out_path / "chr_qc.json", qc)
        return {
            "status": "pass",
            "sample": sample,
            "contig": contig,
            "zarr": str(out_path),
            "shared_read_count": qc.get("shared_read_count", 0),
            "haplotypes": qc.get("haplotypes", {}),
            "read_haplotype_counts": qc.get("read_haplotype_counts", {}),
            "read_phase_set": qc.get("read_phase_set", {}),
            "elapsed_seconds": qc["elapsed_seconds"],
        }
    except BaseException:
        if partial.exists():
            shutil.rmtree(partial)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description="Shared-read HAP1/HAP2 model-input zarr preprocessing from FASTA/chain/BAM TSV")
    parser.add_argument("--record-tsv", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--windows-per-chunk", type=int, default=32)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--schedule", choices=("input", "longest-first"), default="longest-first", help="worker job submission order; longest-first reduces tail latency for multi-sample runs")
    parser.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--window-size", type=int, default=1 << 20)
    parser.add_argument("--offset-unit", type=int, help="deterministic random offset range; default is --window-size, set 1 for fixed-start windows")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--reference-fasta", type=Path, help="Local reference FASTA for CRAM decoding")
    parser.add_argument("--max-records", type=int)
    parser.add_argument("--max-windows", type=int, help="process at most this many eligible windows per HAP; intended for smoke tests")
    parser.add_argument(
        "--validation-level",
        choices=("full", "standard", "fast"),
        default="standard",
        help="QC strictness: full scans read/window strings, standard skips read string scan, fast skips all string scans",
    )
    args = parser.parse_args()
    if args.windows_per_chunk <= 0:
        parser.error("--windows-per-chunk must be positive")
    if args.jobs <= 0:
        parser.error("--jobs must be positive")
    if args.max_windows is not None and args.max_windows <= 0:
        parser.error("--max-windows must be positive")

    config = ExportConfig(
        window_size=args.window_size,
        offset_unit=args.offset_unit,
        seed=args.seed,
        min_mapping_quality=args.min_mapq,
    )
    records = read_record_tsv(args.record_tsv)
    if args.max_records is not None:
        records = records[:args.max_records]
    jobs = [{
        "sample": r["sample"],
        "contig": r["contig"],
        "hap1_fasta": r["hap1_fasta"],
        "hap1_chain": r["hap1_chain"],
        "hap2_fasta": r["hap2_fasta"],
        "hap2_chain": r["hap2_chain"],
        "bam": r["bam"],
        "reference_fasta": str(args.reference_fasta) if args.reference_fasta else "",
        "output_root": str(args.output_root),
        "windows_per_chunk": args.windows_per_chunk,
        "resume": args.resume,
        "overwrite": args.overwrite,
        "max_windows": args.max_windows,
        "validation_level": args.validation_level,
        "config": config.to_dict(),
    } for r in records]
    if args.schedule == "longest-first":
        jobs = sorted(jobs, key=_job_schedule_key)

    manifest_path = args.output_root / "manifest.json"
    manifest = init_manifest(
        record_tsv=args.record_tsv,
        output_root=args.output_root,
        windows_per_chunk=args.windows_per_chunk,
        jobs=args.jobs,
        schedule=args.schedule,
        max_windows=args.max_windows,
        validation_level=args.validation_level,
        config=config.to_dict(),
    )
    write_json(manifest_path, manifest)

    try:
        if args.jobs == 1:
            for job in jobs:
                result = _run_contig(job)
                update_manifest_record(manifest, result)
                write_json(manifest_path, manifest)
                print(json.dumps(result, ensure_ascii=False), flush=True)
        else:
            with ProcessPoolExecutor(max_workers=args.jobs) as pool:
                futures = [pool.submit(_run_contig, job) for job in jobs]
                for future in as_completed(futures):
                    result = future.result()
                    update_manifest_record(manifest, result)
                    write_json(manifest_path, manifest)
                    print(json.dumps(result, ensure_ascii=False), flush=True)
    except BaseException as exc:
        mark_manifest_failed(manifest, exc)
        write_json(manifest_path, manifest)
        raise

    mark_manifest_complete(manifest)
    write_json(manifest_path, manifest)


if __name__ == "__main__":
    main()
