"""Stage 04: uBAM/FASTQ/CRAM/BAM을 canonical coordinate-sorted analysis BAM으로 표준화한다."""
from __future__ import annotations

from collections import defaultdict
import argparse
import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, BinaryIO

import pysam

class PrepareAlignmentError(RuntimeError):
    pass


def run_checked(command: list[str], *, stdout: BinaryIO | int | None = None,
                stderr: BinaryIO | int | None = None) -> subprocess.CompletedProcess[bytes]:
    try:
        return subprocess.run(command, check=True, stdout=stdout, stderr=stderr)
    except subprocess.CalledProcessError as exc:
        raise PrepareAlignmentError(f"명령 실패({exc.returncode}): {' '.join(command)}") from exc


def tool_version(command: list[str]) -> str:
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    except OSError:
        return "not_found"
    lines = (result.stdout or "").splitlines()
    return lines[0] if lines else "unknown"


def bwa_mem2_version(bwa_mem2: str) -> str:
    """Extract BWA-MEM2's semantic version, not its SIMD launcher banner."""
    try:
        result = subprocess.run([bwa_mem2, "version"], stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True)
    except OSError:
        return "bwa-mem2 not_found"
    match = re.search(r"\b\d+(?:\.\d+){1,3}\b", result.stdout or "")
    if match:
        return f"bwa-mem2 {match.group(0)}"
    lines = (result.stdout or "").splitlines()
    return f"bwa-mem2 {lines[-1]}" if lines else "bwa-mem2 unknown"


def require_executable(command: str, label: str) -> None:
    path = Path(command)
    if (path.is_file() and os.access(path, os.X_OK)) or shutil.which(command):
        return
    raise PrepareAlignmentError(f"{label} 실행 파일이 없습니다: {command}")


def normalise_rg(value: str, fallback: str) -> str:
    raw = value.strip() if value and value.strip() != "-" else fallback
    clean = re.sub(r"\s+", "_", raw)
    if not clean:
        raise PrepareAlignmentError("read group ID가 비어 있습니다")
    return clean


def platform_alignment_parameters(platform: str) -> tuple[str, str]:
    """SAM PL과 minimap2 preset을 검증된 플랫폼 이름으로 해석한다."""
    parameters = {
        "hifi": ("PACBIO", "map-hifi"),
        "clr": ("PACBIO", "map-pb"),
        "ont": ("ONT", "map-ont"),
        "illumina": ("ILLUMINA", "sr"),
    }
    try:
        return parameters[platform]
    except KeyError as exc:
        raise PrepareAlignmentError(f"Stage 04 raw long-read 플랫폼 미지원: {platform}") from exc


def index_candidates(path: Path) -> list[Path]:
    if path.suffix == ".cram":
        return [Path(str(path) + ".crai"), path.with_suffix(".crai")]
    return [Path(str(path) + ".bai"), path.with_suffix(".bai")]


def inspect_base_qualities(path: Path, samtools: str, limit: int = 1000) -> dict[str, int]:
    process = subprocess.Popen([samtools, "view", "-F", "4", str(path)],
                               stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    assert process.stdout is not None
    reads = missing = length_mismatch = 0
    quality_values: set[int] = set()
    try:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            reads += 1
            sequence, quality = fields[9], fields[10]
            if sequence == "*" or quality == "*":
                missing += 1
            elif len(sequence) != len(quality):
                length_mismatch += 1
            else:
                quality_values.update(ord(value) - 33 for value in quality)
            if reads >= limit:
                break
    finally:
        process.stdout.close()
        if process.poll() is None:
            process.terminate()
        process.wait()
    return {"reads_inspected": reads, "missing": missing,
            "length_mismatch": length_mismatch,
            "distinct_phred_values": len(quality_values),
            "minimum_phred": min(quality_values) if quality_values else -1,
            "maximum_phred": max(quality_values) if quality_values else -1}


def validate_alignment(path: Path, expected_samples: set[str], samtools: str, *, require_index: bool = True,
                       min_quality_values: int = 2) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size == 0:
        raise PrepareAlignmentError(f"analysis alignment가 없거나 비었습니다: {path}")
    run_checked([samtools, "quickcheck", str(path)], stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE)
    header = subprocess.run([samtools, "view", "-H", str(path)], check=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout
    sq = sum(line.startswith("@SQ\t") for line in header.splitlines())
    if sq == 0:
        raise PrepareAlignmentError(f"좌표 reference @SQ가 없는 alignment입니다: {path}")
    sm_values = set()
    sort_order = None
    for line in header.splitlines():
        fields = line.split("\t")
        if fields and fields[0] == "@HD":
            sort_order = next((x[3:] for x in fields[1:] if x.startswith("SO:")), None)
        if fields and fields[0] == "@RG":
            sm_values.update(x[3:] for x in fields[1:] if x.startswith("SM:"))
    if not expected_samples or not sm_values or not sm_values.issubset(expected_samples):
        raise PrepareAlignmentError(
            f"@RG SM이 허용된 sample alias와 다릅니다: "
            f"expected={sorted(expected_samples)} observed={sorted(sm_values)}"
        )
    if sort_order != "coordinate":
        raise PrepareAlignmentError(f"coordinate sort가 아닙니다: SO={sort_order}")
    if require_index and not any(candidate.is_file() and candidate.stat().st_size > 0
                                 for candidate in index_candidates(path)):
        raise PrepareAlignmentError(f"alignment index가 없습니다: {path}")
    idxstats = subprocess.run([samtools, "idxstats", str(path)], check=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout
    mapped = unmapped = 0
    for line in idxstats.splitlines():
        fields = line.split("\t")
        if len(fields) >= 4:
            mapped += int(fields[2]); unmapped += int(fields[3])
    quality = inspect_base_qualities(path, samtools)
    if quality["reads_inspected"] == 0 or quality["missing"] or quality["length_mismatch"]:
        raise PrepareAlignmentError(f"base quality 누락/길이 오류: {quality}")
    if quality["distinct_phred_values"] < min_quality_values:
        raise PrepareAlignmentError(
            f"base quality가 변별력 없는 placeholder로 의심됩니다: {quality}"
        )
    return {"sq_count": sq, "sort_order": sort_order, "sample_ids": sorted(sm_values),
            "expected_sample_ids": sorted(expected_samples),
            "mapped_reads": mapped, "unmapped_reads": unmapped, "base_quality": quality}


def minimap_command(*, sample: str, platform: str, run_id: str, read_group_id: str,
                  reference: Path, minimap2: str, align_threads: int) -> list[str]:
    """Build the single platform-specific alignment command used by Stage 04."""
    pl, preset = platform_alignment_parameters(platform)
    rg_id = normalise_rg(read_group_id, run_id)
    rg = f"@RG\\tID:{rg_id}\\tSM:{sample}\\tPL:{pl}\\tPU:{run_id}"
    return [minimap2, "-a", "-x", preset, "--MD", "-Y", "-t", str(align_threads),
            "-R", rg, str(reference)]


def bwa_mem2_command(*, sample: str, run_id: str, read_group_id: str, reference: Path,
                     bwa_mem2: str, align_threads: int) -> list[str]:
    """Build the raw Illumina BWA-MEM2 command.

    -Y retains supplementary sequence for the downstream segment contract. -M is
    intentionally omitted because it changes supplementary records into secondary
    records and would blur the later segment-specific QC/filter policy.
    """
    rg_id = normalise_rg(read_group_id, run_id)
    tab = chr(9)
    rg = f"@RG{tab}ID:{rg_id}{tab}SM:{sample}{tab}PL:ILLUMINA{tab}PU:{run_id}"
    return [bwa_mem2, "mem", "-t", str(align_threads), "-Y", "-R", rg, str(reference)]


def bwa_mem2_index_files(prefix: Path) -> list[Path]:
    """Return the complete BWA-MEM2 index file set for one index prefix."""
    return [Path(str(prefix) + suffix) for suffix in
            (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")]


def ensure_bwa_mem2_index(reference: Path, *, bwa_mem2: str, index_root: Path,
                          log: BinaryIO, resource_output: Path | None = None,
                          address_space_bytes: int = 0) -> dict[str, Any]:
    """Build or reuse one locked, reference-fingerprint-specific BWA-MEM2 index.

    The source FASTA directory is never mutated. The completion manifest is written
    only after every component has been moved under the shared cache prefix.
    """
    require_executable(bwa_mem2, "BWA-MEM2")
    fingerprint = reference_fingerprint(reference)
    key = hashlib.sha256(json.dumps(fingerprint, sort_keys=True).encode()).hexdigest()[:24]
    cache_dir = index_root / key
    cache_dir.mkdir(parents=True, exist_ok=True)
    prefix = cache_dir / "reference"
    complete = cache_dir / "complete.json"
    lock_path = cache_dir / ".build.lock"
    with lock_path.open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            if complete.is_file() and all(path.is_file() and path.stat().st_size > 0
                                          for path in bwa_mem2_index_files(prefix)):
                try:
                    metadata = json.loads(complete.read_text())
                except (OSError, json.JSONDecodeError):
                    metadata = {}
                if metadata.get("reference") == fingerprint:
                    return {**metadata, "status": "reused"}

            build_dir = Path(tempfile.mkdtemp(prefix=".build.", dir=cache_dir))
            temporary_prefix = build_dir / "reference"
            try:
                command = [bwa_mem2, "index", "-p", str(temporary_prefix), str(reference)]
                if resource_output is not None:
                    command = [sys.executable, str(Path(__file__).resolve().parents[2] / 'scripts/resource_runner.py'),
                               '--output', str(resource_output), '--address-space-bytes', str(address_space_bytes),
                               '--', *command]
                run_checked(command,
                            stdout=subprocess.DEVNULL, stderr=log)
                missing = [path for path in bwa_mem2_index_files(temporary_prefix)
                           if not path.is_file() or path.stat().st_size == 0]
                if missing:
                    raise PrepareAlignmentError(f"BWA-MEM2 index 생성 불완전: {missing}")
                for source, destination in zip(bwa_mem2_index_files(temporary_prefix),
                                               bwa_mem2_index_files(prefix)):
                    os.replace(source, destination)
                metadata = {"reference": fingerprint, "index_prefix": str(prefix.resolve()),
                            "files": [str(path.resolve()) for path in bwa_mem2_index_files(prefix)]}
                write_json_atomic(complete, metadata)
                return {**metadata, "status": "created"}
            finally:
                shutil.rmtree(build_dir, ignore_errors=True)
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def sort_aligner_output(aligner: subprocess.Popen[bytes], output: Path, *, samtools: str,
                        sort_threads: int, sort_memory: str, log: BinaryIO) -> int:
    """Coordinate-sort a minimap2 SAM stream into one temporary BAM."""
    assert aligner.stdout is not None
    sorter = subprocess.Popen(
        [samtools, "sort", "-@", str(sort_threads), "-m", sort_memory, "-T",
         str(output) + ".sort", "-o", str(output), "-"],
        stdin=aligner.stdout, stdout=subprocess.DEVNULL, stderr=log,
    )
    aligner.stdout.close()
    sort_rc = sorter.wait()
    align_rc = aligner.wait()
    if (align_rc, sort_rc) != (0, 0):
        raise PrepareAlignmentError(f"aligner/sort 실패: aligner={align_rc} sort={sort_rc}")
    return align_rc


def align_long_read_row(row: dict[str, str], output: Path, *, sample: str, platform: str,
                        reference: Path, minimap2: str, samtools: str, align_threads: int,
                        sort_threads: int, sort_memory: str, extract_threads: int,
                        log: BinaryIO) -> None:
    """Align one raw ONT/HiFi/CLR uBAM or FASTQ into a sorted BAM."""
    fmt = row["format"].lower()
    source = Path(row["path"])
    command = minimap_command(sample=sample, platform=platform, run_id=row["run_id"],
                              read_group_id=row.get("read_group_id", ""), reference=reference,
                              minimap2=minimap2, align_threads=align_threads)
    extractor: subprocess.Popen[bytes] | None = None
    if fmt == "ubam":
        extractor = subprocess.Popen([samtools, "fastq", "-@", str(extract_threads), "-n", str(source)],
                                     stdout=subprocess.PIPE, stderr=log)
        assert extractor.stdout is not None
        aligner = subprocess.Popen(command + ["-"], stdin=extractor.stdout, stdout=subprocess.PIPE, stderr=log)
        extractor.stdout.close()
    elif fmt in {"fastq", "fastq_gz"}:
        aligner = subprocess.Popen(command + [str(source)], stdout=subprocess.PIPE, stderr=log)
    else:
        raise PrepareAlignmentError(f"Stage 04 raw long-read 형식 미지원: {fmt}")
    try:
        sort_aligner_output(aligner, output, samtools=samtools, sort_threads=sort_threads,
                            sort_memory=sort_memory, log=log)
    finally:
        if extractor is not None:
            extract_rc = extractor.wait()
            if extract_rc != 0:
                raise PrepareAlignmentError(f"uBAM FASTQ 추출 실패: run={row['run_id']} rc={extract_rc}")


def align_illumina_queries(query_files: list[Path], output: Path, *, sample: str, run_id: str,
                            read_group_id: str, reference: Path, minimap2: str, bwa_mem2: str,
                            bwa_mem2_reference: Path | None, illumina_aligner: str, samtools: str,
                            align_threads: int, sort_threads: int, sort_memory: str,
                            log: BinaryIO) -> None:
    """Align one paired R1/R2 set or one singleton stream with the explicit backend."""
    if not query_files or len(query_files) > 2:
        raise PrepareAlignmentError(f"Illumina query file 수 오류: run={run_id} files={query_files}")
    if illumina_aligner == "bwa-mem2":
        command = bwa_mem2_command(sample=sample, run_id=run_id, read_group_id=read_group_id,
                                   reference=bwa_mem2_reference or reference, bwa_mem2=bwa_mem2,
                                   align_threads=align_threads)
    elif illumina_aligner == "minimap2":
        command = minimap_command(sample=sample, platform="illumina", run_id=run_id,
                                  read_group_id=read_group_id, reference=reference,
                                  minimap2=minimap2, align_threads=align_threads)
    else:
        raise PrepareAlignmentError(f"지원하지 않는 Illumina aligner: {illumina_aligner}")
    aligner = subprocess.Popen(command + [str(path) for path in query_files],
                               stdout=subprocess.PIPE, stderr=log)
    sort_aligner_output(aligner, output, samtools=samtools, sort_threads=sort_threads,
                        sort_memory=sort_memory, log=log)


def raw_illumina_bams(files: list[dict[str, str]], temp_dir: Path, *, sample: str, reference: Path,
                      minimap2: str, bwa_mem2: str, illumina_aligner: str, samtools: str,
                      bwa_mem2_reference: Path | None, align_threads: int, sort_threads: int,
                      sort_memory: str, extract_threads: int, log: BinaryIO) -> list[Path]:
    """Create sorted BAMs for every raw Illumina run, including explicit singletons."""
    if illumina_aligner == "bwa-mem2":
        require_executable(bwa_mem2, "BWA-MEM2")
        if bwa_mem2_reference is None:
            raise PrepareAlignmentError("BWA-MEM2 index prefix가 준비되지 않았습니다")
    elif illumina_aligner == "minimap2":
        require_executable(minimap2, "minimap2")
    else:
        raise PrepareAlignmentError(f"지원하지 않는 Illumina aligner: {illumina_aligner}")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in files:
        grouped[row["run_id"]].append(row)
    outputs: list[Path] = []
    ordinal = 0
    for run_id, rows in sorted(grouped.items()):
        by_role: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            by_role[row["file_role"].lower()].append(row)
        rg_values = {normalise_rg(row.get("read_group_id", ""), run_id) for row in rows}
        if len(rg_values) != 1:
            raise PrepareAlignmentError(f"같은 Illumina run의 read_group_id가 다릅니다: run={run_id}")
        rg_id = next(iter(rg_values))
        paired = by_role.get("read1", []) + by_role.get("read2", [])
        if paired:
            if len(by_role.get("read1", [])) != 1 or len(by_role.get("read2", [])) != 1:
                raise PrepareAlignmentError(f"Illumina R1/R2 입력 오류: run={run_id}")
            ordinal += 1
            output = temp_dir / f"run_{ordinal:04d}.paired.bam"
            align_illumina_queries([Path(by_role["read1"][0]["path"]), Path(by_role["read2"][0]["path"])],
                                    output, sample=sample, run_id=run_id, read_group_id=rg_id,
                                    reference=reference, minimap2=minimap2, bwa_mem2=bwa_mem2,
                                    bwa_mem2_reference=bwa_mem2_reference,
                                    illumina_aligner=illumina_aligner, samtools=samtools,
                                    align_threads=align_threads, sort_threads=sort_threads,
                                    sort_memory=sort_memory, log=log)
            outputs.append(output)
        for row in by_role.get("singletons", []):
            ordinal += 1
            output = temp_dir / f"run_{ordinal:04d}.singleton.bam"
            align_illumina_queries([Path(row["path"])], output, sample=sample, run_id=run_id,
                                    read_group_id=rg_id, reference=reference, minimap2=minimap2,
                                    bwa_mem2=bwa_mem2, bwa_mem2_reference=bwa_mem2_reference,
                                    illumina_aligner=illumina_aligner, samtools=samtools,
                                    align_threads=align_threads,
                                    sort_threads=sort_threads, sort_memory=sort_memory, log=log)
            outputs.append(output)
        for row in by_role.get("alignment", []):
            if row["format"].lower() != "ubam":
                raise PrepareAlignmentError(f"raw Illumina alignment은 uBAM만 허용: run={run_id}")
            prefix = temp_dir / f"run_{run_id}.ubam"
            r1, r2, single = (Path(str(prefix) + suffix) for suffix in (".R1.fastq", ".R2.fastq", ".single.fastq"))
            run_checked([samtools, "fastq", "-@", str(extract_threads), "-n", "-1", str(r1), "-2", str(r2),
                         "-0", "/dev/null", "-s", str(single), str(row["path"])], stderr=log)
            for label, queries in (("paired", [r1, r2]), ("singleton", [single])):
                if not all(path.is_file() and path.stat().st_size for path in queries):
                    continue
                ordinal += 1
                output = temp_dir / f"run_{ordinal:04d}.ubam_{label}.bam"
                align_illumina_queries(queries, output, sample=sample, run_id=run_id,
                                        read_group_id=rg_id, reference=reference, minimap2=minimap2,
                                        bwa_mem2=bwa_mem2, bwa_mem2_reference=bwa_mem2_reference,
                                        illumina_aligner=illumina_aligner, samtools=samtools,
                                        align_threads=align_threads,
                                        sort_threads=sort_threads, sort_memory=sort_memory, log=log)
                outputs.append(output)
    if not outputs:
        raise PrepareAlignmentError("raw Illumina FASTQ/uBAM에서 정렬할 read가 없습니다")
    return outputs


def reference_fingerprint(reference: Path) -> dict[str, Any]:
    stat = reference.stat()
    fai = Path(str(reference) + ".fai")
    return {"path": str(reference.resolve()), "size_bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns, "fai_path": str(fai.resolve()) if fai.is_file() else None,
            "fai_size_bytes": fai.stat().st_size if fai.is_file() else None}


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def require_reference(path: Path | None, label: str) -> Path:
    if path is None or not path.is_file() or not Path(str(path) + ".fai").is_file():
        raise PrepareAlignmentError(f"{label} 또는 .fai가 없습니다: {path}")
    return path


def source_rows(contract: dict[str, Any]) -> list[dict[str, str]]:
    """Normalize legacy and file-manifest aligned/raw inputs into manifest-like rows."""
    source = contract["source"]
    if contract["input_mode"] == "legacy_alignment":
        path = str(source.get("alignment") or "")
        if not path:
            raise PrepareAlignmentError("legacy alignment 경로가 비어 있습니다")
        return [{"run_id": "legacy", "file_role": "alignment",
                 "format": Path(path).suffix.lower().lstrip("."), "alignment_state": "aligned",
                 "path": path, "source_reference_fasta": "", "read_group_id": "legacy"}]
    return [dict(row) for row in source.get("files", [])]


def source_identity(row: dict[str, str]) -> dict[str, Any]:
    path = Path(row["path"]).resolve()
    if not path.is_file():
        raise PrepareAlignmentError(f"source input이 없습니다: {path}")
    stat = path.stat()
    return {"path": str(path), "format": row["format"].lower(), "run_id": row["run_id"],
            "file_role": row["file_role"], "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def source_reference(row: dict[str, str], fallback: Path | None) -> Path | None:
    raw = row.get("source_reference_fasta") or ""
    return Path(raw) if raw else fallback


def _tag_record_stats_python(path: Path, reference: Path | None = None) -> dict[str, int]:
    """Full-record counters used to prove conversion/marking did not lose records or tags."""
    kwargs: dict[str, str] = {}
    if path.suffix.lower() == ".cram":
        if reference is None:
            raise PrepareAlignmentError(f"CRAM decode reference가 없습니다: {path}")
        kwargs["reference_filename"] = str(reference)
    counts = {
        "records": 0, "mapped_records": 0, "primary_records": 0,
        "secondary_records": 0, "supplementary_records": 0,
        "duplicate_records": 0, "primary_duplicate_records": 0,
        "hp_tagged_records": 0, "ps_tagged_records": 0,
    }
    with pysam.AlignmentFile(str(path), "rb", **kwargs) as alignment:
        for read in alignment.fetch(until_eof=True):
            counts["records"] += 1
            if not read.is_unmapped:
                counts["mapped_records"] += 1
            if read.is_secondary:
                counts["secondary_records"] += 1
            elif read.is_supplementary:
                counts["supplementary_records"] += 1
            else:
                counts["primary_records"] += 1
            if read.is_duplicate:
                counts["duplicate_records"] += 1
                if not read.is_secondary and not read.is_supplementary:
                    counts["primary_duplicate_records"] += 1
            if read.has_tag("HP"):
                counts["hp_tagged_records"] += 1
            if read.has_tag("PS"):
                counts["ps_tagged_records"] += 1
    return counts


def tag_record_stats(path: Path, reference: Path | None = None, *, samtools: str | None = None,
                     details: dict | None = None) -> dict[str, int]:
    from integrate_wgs_rna.prepare_phases import mark_read_only_prepare
    mark_read_only_prepare()
    backend = os.environ.get("MFM_RECORD_STATS_BACKEND", "auto")
    if backend not in {"auto", "python", "native"}:
        raise PrepareAlignmentError("MFM_RECORD_STATS_BACKEND는 auto/python/native여야 함")
    if backend == "python" or (backend == "auto" and path.stat().st_size < 64 * 1024**2):
        return _tag_record_stats_python(path, reference)
    from integrate_wgs_rna.alignment_record_stats import native_record_stats
    return native_record_stats(path, reference, samtools=samtools, details=details)


def illumina_duplicate_provenance(row: dict[str, str], fallback_reference: Path | None) -> dict[str, Any]:
    """Accept only Picard/GATK MarkDuplicates provenance that did not remove records."""
    path = Path(row["path"])
    reference = source_reference(row, fallback_reference)
    if path.suffix.lower() == ".cram":
        reference = require_reference(reference, "CRAM source decode reference")
    kwargs: dict[str, str] = {}
    if path.suffix.lower() == ".cram":
        kwargs["reference_filename"] = str(reference)
    with pysam.AlignmentFile(str(path), "rb", **kwargs) as alignment:
        programs = alignment.header.to_dict().get("PG", [])
    matches: list[dict[str, str]] = []
    removed = False
    for program in programs:
        joined = " ".join(str(program.get(key, "")) for key in ("ID", "PN", "CL")).lower()
        if "markduplicates" not in joined:
            continue
        command = str(program.get("CL", "")).lower()
        explicitly_removed = bool(re.search(
            r"\bremove_(?:sequencing_)?duplicates\s*(?:=|\s)\s*true\b", command))
        removed = removed or explicitly_removed
        match = {key: str(program[key]) for key in ("ID", "PN", "VN") if key in program}
        match["records_removed"] = str(explicitly_removed).lower()
        matches.append(match)
    return {
        "path": str(path.resolve()),
        "trusted_picard_marking": bool(matches) and not removed,
        "matched_programs": matches,
    }


def resolve_duplicate_plan(*, platform: str, states: list[str], requested_policy: str,
                           provenance: list[dict[str, Any]]) -> dict[str, Any]:
    """Choose one canonical-coordinate duplicate decision before personal-HAP alignment."""
    if requested_policy not in {"auto", "mark", "preserve"}:
        raise PrepareAlignmentError(
            f"Illumina duplicate policy는 auto|mark|preserve여야 합니다: {requested_policy}")
    if platform != "illumina":
        return {"requested_policy": requested_policy, "action": "not_applicable",
                "reason": "long-read duplicate marking은 library/tool provenance 의존적이므로 자동 적용하지 않음"}
    trusted = bool(provenance) and all(item.get("trusted_picard_marking") for item in provenance)
    if requested_policy == "preserve":
        if states != ["aligned"] or not trusted:
            raise PrepareAlignmentError(
                "duplicate preserve는 모든 aligned Illumina source의 Picard MarkDuplicates 이력이 필요합니다")
        action = "preserved_source_picard"
    elif requested_policy == "mark":
        action = "picard_marked"
    else:
        action = "preserved_source_picard" if states == ["aligned"] and trusted else "picard_marked"
    return {"requested_policy": requested_policy, "action": action,
            "source_provenance": provenance, "remove_duplicates": False}


def mark_illumina_duplicates(input_bam: Path, marked_bam: Path, metrics_path: Path, *,
                             picard: str, heap_gb: int, temp_dir: Path, log: BinaryIO) -> None:
    """Picard-mark Illumina duplicates while retaining every alignment record."""
    if heap_gb < 1:
        raise PrepareAlignmentError(f"Picard heap GB는 양의 정수여야 합니다: {heap_gb}")
    require_executable(picard, "Picard")
    picard_tmp = temp_dir / "picard_tmp"
    picard_tmp.mkdir(parents=True, exist_ok=True)
    run_checked([
        picard, f"-Xmx{heap_gb}g", f"-Djava.io.tmpdir={picard_tmp}", "MarkDuplicates",
        f"I={input_bam}", f"O={marked_bam}", f"M={metrics_path}",
        "ASSUME_SORT_ORDER=coordinate", "VALIDATION_STRINGENCY=LENIENT",
        "REMOVE_DUPLICATES=false", f"TMP_DIR={picard_tmp}", "CREATE_INDEX=false",
    ], stdout=log, stderr=log)


def merge_sorted_bams(bams: list[Path], merged: Path, *, samtools: str, sort_threads: int,
                      log: BinaryIO) -> None:
    if not bams:
        raise PrepareAlignmentError("병합할 BAM이 없습니다")
    if len(bams) == 1:
        os.replace(bams[0], merged)
        return
    run_checked([samtools, "merge", "-@", str(sort_threads), "-f", str(merged), *map(str, bams)],
                stdout=subprocess.DEVNULL, stderr=log)


def transcode_aligned_source(row: dict[str, str], output: Path, *, samtools: str, sort_threads: int,
                             sort_memory: str, fallback_reference: Path | None, log: BinaryIO) -> dict[str, int]:
    """Decode one aligned BAM/CRAM and coordinate-sort it into canonical BAM form."""
    source = Path(row["path"])
    fmt = row["format"].lower()
    if fmt not in {"bam", "cram"}:
        raise PrepareAlignmentError(f"aligned source 형식 오류: {fmt}")
    reference = source_reference(row, fallback_reference)
    if fmt == "cram":
        reference = require_reference(reference, "CRAM source decode reference")
    run_checked([samtools, "quickcheck", str(source)], stdout=subprocess.DEVNULL, stderr=log)
    details = {}
    before = tag_record_stats(source, reference, samtools=samtools, details=details)
    from integrate_wgs_rna.prepare_phases import reserve_prepare_write
    reserve_prepare_write()
    from integrate_wgs_rna.cpu_executor import WorkerBudget
    from integrate_wgs_rna.input_parallel import cpu_budget
    # 전수 집계가 반환한 CPU를 실제 view/sort 시작 전에 다시 예약함.
    total = cpu_budget()
    budget = WorkerBudget(total, total, True)
    budget.base_memory = 64 * 1024**3
    from integrate_wgs_rna.gpu_bam_io import enabled as gpu_io_enabled
    gpu_io = gpu_io_enabled()
    minimum = 6 if gpu_io else 1
    if total < minimum:
        raise PrepareAlignmentError("nvcomp BAM 준비에는 최소 6CPU 예약이 필요함")
    admitted = budget.capacity(total)
    while admitted < minimum:
        import time
        time.sleep(.2)
        admitted = budget.capacity(total)
    io_cpu = 3 if gpu_io else 0
    if (not gpu_io and details.get('coordinate_sorted') is True
            and details.get('header_coordinate') is True):
        from integrate_wgs_rna.alignment_record_stats import fingerprint
        if fingerprint(source) != details['input_fingerprint']:
            raise PrepareAlignmentError('정렬 순서 확인 뒤 입력이 변경됨')
        from integrate_wgs_rna.prepare_memory import reserve_transcode
        reserve_transcode(admitted, 0, sort_memory)
        command = [samtools, 'view', '-@', str(max(0, admitted-2)), '-O', 'bam,level=1']
        if fmt == 'cram':command += ['-T', str(reference)]
        print('[ALIGNED-TRANSCODE-PLAN] '+json.dumps(dict(cpus=admitted,
              coordinate_order_verified=True, sort_skipped=True, io_backend='cpu')), flush=True)
        run_checked([*command, '-o', str(output), str(source)],
                    stdout=subprocess.DEVNULL, stderr=log)
        if fingerprint(source) != details['input_fingerprint']:
            raise PrepareAlignmentError('빠른 변환 중 입력이 변경됨')
        return before
    sort_threads = min(sort_threads, max(0, (admitted - 2 - io_cpu) // 2))
    view_threads = max(0, admitted - sort_threads - 2 - io_cpu)
    view_command = [samtools, "view", "-@", str(view_threads), "-u"]
    if fmt == "cram":
        view_command.extend(["-T", str(reference)])
    view_command.append(str(source))
    print("[ALIGNED-TRANSCODE-PLAN] " + json.dumps({"cpus": admitted,
          "view_additional_threads": view_threads, "sort_additional_threads": sort_threads,
          "single_cpu_sequential": admitted < 2, "io_backend": "nvcomp" if gpu_io else "cpu"}), flush=True)
    if gpu_io:
        from integrate_wgs_rna.gpu_bam_io import sort_stream
        sort_stream(view_command, output, samtools=samtools, sort_threads=sort_threads,
                    sort_memory=sort_memory, log=log)
        return before
    from integrate_wgs_rna.prepare_memory import reserve_transcode
    reserve_transcode(admitted, sort_threads, sort_memory)
    if admitted < 2:
        decoded = output.with_name(output.name + ".decoded.part.bam")
        try:
            decode_command = ["-b" if item == "-u" else item for item in view_command[:-1]]
            run_checked([*decode_command, "-o", str(decoded), str(source)],
                        stdout=subprocess.DEVNULL, stderr=log)
            run_checked([samtools, "sort", "-@", "0", "-m", sort_memory,
                         "-T", str(output) + ".sort", "-o", str(output), str(decoded)],
                        stdout=subprocess.DEVNULL, stderr=log)
        finally:
            decoded.unlink(missing_ok=True)
        return before
    viewer = subprocess.Popen(view_command, stdout=subprocess.PIPE, stderr=log)
    assert viewer.stdout is not None
    sorter = subprocess.Popen(
        [samtools, "sort", "-@", str(sort_threads), "-m", sort_memory, "-T",
         str(output) + ".sort", "-o", str(output), "-"],
        stdin=viewer.stdout, stdout=subprocess.DEVNULL, stderr=log,
    )
    viewer.stdout.close()
    sort_rc = sorter.wait()
    view_rc = viewer.wait()
    if (view_rc, sort_rc) != (0, 0):
        raise PrepareAlignmentError(f"aligned→BAM 변환 실패: view={view_rc} sort={sort_rc}")
    return before


def finalise_canonical_bam(merged: Path, output: Path, *, expected_samples: set[str], samtools: str, sort_threads: int,
                           min_quality_values: int) -> tuple[dict[str, Any], Path]:
    from integrate_wgs_rna.prepare_phases import mark_read_only_prepare
    mark_read_only_prepare(index_threads=sort_threads)
    part = output.with_name(output.name + f".part.{os.getpid()}")
    part_index = Path(str(part) + ".bai")
    os.replace(merged, part)
    run_checked([samtools, "index", "-@", str(sort_threads), str(part), str(part_index)],
                stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    metrics = validate_alignment(part, expected_samples, samtools, min_quality_values=min_quality_values)
    final_index = Path(str(output) + ".bai")
    os.replace(part, output)
    os.replace(part_index, final_index)
    return metrics, final_index


def alignment_policy(*, states: list[str], platform: str, illumina_aligner: str) -> dict[str, str]:
    """Describe the method that made the canonical BAM; use it as a reuse guard."""
    if states == ["aligned"]:
        return {"mode": "container_standardization", "aligner": "none"}
    if platform == "illumina":
        if illumina_aligner == "bwa-mem2":
            return {"mode": "raw_illumina", "aligner": "bwa-mem2", "preset": "mem",
                    "options": "-Y; no -M"}
        if illumina_aligner == "minimap2":
            return {"mode": "raw_illumina", "aligner": "minimap2", "preset": "sr",
                    "options": "--MD -Y"}
        raise PrepareAlignmentError(f"지원하지 않는 Illumina aligner: {illumina_aligner}")
    if platform in {"hifi", "clr", "ont"}:
        return {"mode": "raw_long_read", "aligner": "minimap2",
                "preset": platform_alignment_parameters(platform)[1], "options": "--MD -Y"}
    raise PrepareAlignmentError(f"Stage 04 플랫폼 미지원: {platform}")


def reuse_matches_contract(qc_json: Path, inputs: list[dict[str, Any]], reference: Path,
                           policy: dict[str, str]) -> bool:
    if not qc_json.is_file():
        return False
    try:
        payload = json.loads(qc_json.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    return (payload.get("canonical_format") == "bam"
            and payload.get("source_inputs") == inputs
            and payload.get("normalization_reference") == reference_fingerprint(reference)
            and payload.get("alignment_policy") == policy)


def prepare(contract: dict[str, Any], *, output: Path, qc_json: Path, minimap2: str,
            samtools: str, align_threads: int, sort_threads: int, sort_memory: str,
            extract_threads: int, overwrite: bool, keep_run_bams: bool, min_quality_values: int = 2,
            bwa_mem2: str = "bwa-mem2", illumina_aligner: str = "bwa-mem2",
            bwa_mem2_index_root: Path | None = None, picard: str = "picard",
            illumina_duplicate_policy: str = "auto", markdup_heap_gb: int = 16) -> dict[str, Any]:
    """Create exactly one validated canonical BAM, regardless of the source container."""
    sample = contract["sample_id"]
    source_sample_id = str(contract.get("source_sample_id") or sample).strip()
    if not source_sample_id:
        raise PrepareAlignmentError("source_sample_id가 비어 있습니다")
    expected_sample_ids = {sample, source_sample_id}
    platform = contract["platform"]
    library_provenance = contract.get("library_provenance", {
        "wgs": {"amplification": "unspecified", "protocol": None, "evidence": None},
        "rna": {"amplification": "unspecified", "protocol": None,
                "pcr_cycles": None, "evidence": None},
    })
    illumina_aligner = illumina_aligner.lower()
    states = contract["source"].get("alignment_states", [])
    rows = source_rows(contract)
    if output.suffix.lower() != ".bam":
        raise PrepareAlignmentError(f"Stage 04 canonical output은 BAM이어야 합니다: {output}")
    if not rows:
        raise PrepareAlignmentError("raw/aligned source 행이 없습니다")
    if states not in (["aligned"], ["unaligned"]):
        raise PrepareAlignmentError(f"지원하지 않는 alignment state: {states}")
    policy = alignment_policy(states=states, platform=platform, illumina_aligner=illumina_aligner)
    inputs = [source_identity(row) for row in rows]
    refs = contract.get("references", {})
    source_decode = Path(refs["source_decode"]) if refs.get("source_decode") else None
    analysis_reference = Path(refs["analysis"]) if refs.get("analysis") else None
    if states == ["unaligned"]:
        reference = require_reference(analysis_reference, "analysis reference")
    else:
        reference = require_reference(source_decode or analysis_reference, "source/decode reference")
    source_duplicate_provenance = []
    if platform == "illumina" and states == ["aligned"]:
        source_duplicate_provenance = [
            illumina_duplicate_provenance(row, source_decode) for row in rows
        ]
    duplicate_plan = resolve_duplicate_plan(
        platform=platform, states=states, requested_policy=illumina_duplicate_policy,
        provenance=source_duplicate_provenance,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not overwrite:
        if not reuse_matches_contract(qc_json, inputs, reference, policy):
            raise PrepareAlignmentError(
                f"기존 canonical BAM의 source/reference provenance가 다릅니다: {output}; --overwrite 필요"
            )
        metrics = validate_alignment(output, expected_sample_ids, samtools, min_quality_values=min_quality_values)
        # Raw/unprovenanced Illumina의 과거 canonical BAM은 duplicate marking 없이 재사용하지 않는다.
        prior_qc = json.loads(qc_json.read_text())
        prior_duplicate = prior_qc.get("duplicate_marking", {})
        if (duplicate_plan["action"] == "picard_marked"
                and prior_duplicate.get("action") != "picard_marked"):
            raise PrepareAlignmentError(
                f"기존 Illumina canonical BAM에 Picard duplicate marking QC가 없습니다: {output}; --overwrite 필요"
            )
        output_stats = tag_record_stats(output, samtools=samtools)
        retained_source_stats = prior_qc.get("source_record_tag_counts")
        if retained_source_stats:
            expected = {key: sum(item[key] for item in retained_source_stats) for key in output_stats}
            compare = set(output_stats)
            if prior_duplicate.get("action") == "picard_marked":
                compare -= {"duplicate_records", "primary_duplicate_records"}
            if any(expected[key] != output_stats[key] for key in compare):
                raise PrepareAlignmentError("재사용 BAM의 원본 record/tag 보존 검사 실패")
        duplicate_marking = {
            **duplicate_plan, "status": "pass",
            "before": prior_duplicate.get("before"),
            "after": output_stats,
            "metrics_path": prior_duplicate.get("metrics_path"),
        }
        # policy match를 확인한 기존 QC의 tool/index provenance는 canonical BAM 재사용에도 보존한다.
        payload = {
            "schema": "mfm_stage04_alignment_qc_v4", "sample_id": sample,
            "source_sample_id": source_sample_id,
            "status": "reused_canonical_bam", "canonical_format": "bam",
            "platform": platform, "preset": prior_qc.get("preset"),
            "analysis_alignment": str(output.resolve()), "analysis_index": str(Path(str(output) + ".bai").resolve()),
            "normalization_reference": reference_fingerprint(reference), "source_inputs": inputs,
            "alignment_policy": policy, "bwa_mem2_index": prior_qc.get("bwa_mem2_index"),
            "library_provenance": library_provenance,
            "duplicate_marking": duplicate_marking,
            "tools": prior_qc.get("tools"),
            "output_record_tag_counts": output_stats, "metrics": metrics,
            "source_record_tag_counts": retained_source_stats,
            "source_file_count": prior_qc.get("source_file_count"),
            "source_runs": prior_qc.get("source_runs"), "source_io_execution": prior_qc.get("source_io_execution", []),
        }
        write_json_atomic(qc_json, payload)
        return payload

    temp_dir = Path(tempfile.mkdtemp(prefix=f".stage04.{sample}.", dir=output.parent))
    bwa_mem2_index: dict[str, Any] | None = None
    log_path = output.with_suffix(output.suffix + ".stage04.log")
    try:
        with log_path.open("ab") as log:
            run_bams: list[Path] = []
            source_stats: list[dict[str, int]] = []
            source_io_execution = []
            if states == ["aligned"]:
                for index, row in enumerate(rows, 1):
                    run_bam = temp_dir / f"aligned_{index:04d}.bam"
                    source_stats.append(transcode_aligned_source(
                        row, run_bam, samtools=samtools, sort_threads=sort_threads, sort_memory=sort_memory,
                        fallback_reference=source_decode, log=log,
                    ))
                    io_report = run_bam.with_name(run_bam.name+'.gpu_io.json')
                    if io_report.exists():source_io_execution.append(json.loads(io_report.read_text()))
                    run_bams.append(run_bam)
                status = "converted_aligned"
                preset = "container_standardization"
            elif platform == "illumina":
                from integrate_wgs_rna.prepare_phases import reserve_prepare_write
                reserve_prepare_write()
                bwa_mem2_reference = None
                if illumina_aligner == "bwa-mem2":
                    index_root = bwa_mem2_index_root or (output.parent / "_reference_cache" / "bwa_mem2")
                    bwa_mem2_index = ensure_bwa_mem2_index(
                        reference, bwa_mem2=bwa_mem2, index_root=index_root, log=log)
                    bwa_mem2_reference = Path(bwa_mem2_index["index_prefix"])
                run_bams = raw_illumina_bams(
                    rows, temp_dir, sample=sample, reference=reference, minimap2=minimap2,
                    bwa_mem2=bwa_mem2, illumina_aligner=illumina_aligner, samtools=samtools,
                    bwa_mem2_reference=bwa_mem2_reference,
                    align_threads=align_threads, sort_threads=sort_threads, sort_memory=sort_memory,
                    extract_threads=extract_threads, log=log,
                )
                status = "created_raw_illumina"
                preset = policy["preset"]
            elif platform in {"hifi", "clr", "ont"}:
                from integrate_wgs_rna.prepare_phases import reserve_prepare_write
                reserve_prepare_write()
                for index, row in enumerate(rows, 1):
                    run_bam = temp_dir / f"raw_{index:04d}.bam"
                    align_long_read_row(row, run_bam, sample=sample, platform=platform, reference=reference,
                                        minimap2=minimap2, samtools=samtools, align_threads=align_threads,
                                        sort_threads=sort_threads, sort_memory=sort_memory,
                                        extract_threads=extract_threads, log=log)
                    run_bams.append(run_bam)
                status = "created_raw_longread"
                preset = platform_alignment_parameters(platform)[1]
            else:
                raise PrepareAlignmentError(f"Stage 04 플랫폼 미지원: {platform}")

            merged = temp_dir / f"{sample}.merged.bam"
            from integrate_wgs_rna.prepare_phases import reserve_prepare_write
            reserve_prepare_write()
            merge_sorted_bams(run_bams, merged, samtools=samtools, sort_threads=sort_threads, log=log)
            pre_mark_stats = tag_record_stats(merged, samtools=samtools)
            post_mark_stats = pre_mark_stats
            markdup_metrics_temp: Path | None = None
            markdup_metrics_final: Path | None = None
            if duplicate_plan["action"] == "picard_marked":
                reserve_prepare_write()
                marked = temp_dir / f"{sample}.markduplicates.bam"
                markdup_metrics_temp = temp_dir / f"{sample}.markduplicates.metrics.txt"
                mark_illumina_duplicates(
                    merged, marked, markdup_metrics_temp, picard=picard, heap_gb=markdup_heap_gb,
                    temp_dir=temp_dir, log=log,
                )
                post_mark_stats = tag_record_stats(marked, samtools=samtools)
                invariant_keys = set(pre_mark_stats) - {"duplicate_records", "primary_duplicate_records"}
                if any(pre_mark_stats[key] != post_mark_stats[key] for key in invariant_keys):
                    raise PrepareAlignmentError(
                        f"Picard MarkDuplicates record/tag 보존 실패: before={pre_mark_stats} after={post_mark_stats}")
                os.replace(marked, merged)
            metrics, final_index = finalise_canonical_bam(
                merged, output, expected_samples=expected_sample_ids, samtools=samtools, sort_threads=sort_threads,
                min_quality_values=min_quality_values,
            )
            # finalise_canonical_bam은 검증·index 뒤 같은 BAM을 원자적으로 rename하므로
            # 대용량 파일을 다시 전수 scan하지 않고 직전 검증된 통계를 사용한다.
            output_stats = post_mark_stats
            expected = {key: sum(item[key] for item in source_stats) for key in output_stats} if source_stats else {}
            compare_keys = set(output_stats)
            if duplicate_plan["action"] == "picard_marked":
                compare_keys -= {"duplicate_records", "primary_duplicate_records"}
            if expected and any(expected[key] != output_stats[key] for key in compare_keys):
                output.unlink(missing_ok=True)
                Path(str(output) + ".bai").unlink(missing_ok=True)
                raise PrepareAlignmentError(
                    f"aligned→BAM record/tag 보존 실패: source={expected} output={output_stats}")
            if markdup_metrics_temp is not None:
                markdup_metrics_final = output.with_suffix(output.suffix + ".markduplicates.metrics.txt")
                os.replace(markdup_metrics_temp, markdup_metrics_final)
            duplicate_marking = {
                **duplicate_plan, "status": "pass", "before": pre_mark_stats,
                "after": output_stats,
                "metrics_path": str(markdup_metrics_final.resolve()) if markdup_metrics_final else None,
            }

        tools = {"samtools": tool_version([samtools, "--version"])}
        if policy["aligner"] == "minimap2":
            tools["minimap2"] = tool_version([minimap2, "--version"])
        elif policy["aligner"] == "bwa-mem2":
            tools["bwa_mem2"] = bwa_mem2_version(bwa_mem2)
        if duplicate_plan["action"] == "picard_marked":
            tools["picard"] = tool_version([picard, "MarkDuplicates", "--version"])
        payload = {
            "schema": "mfm_stage04_alignment_qc_v4", "sample_id": sample,
            "source_sample_id": source_sample_id, "status": status,
            "canonical_format": "bam", "platform": platform, "preset": preset,
            "analysis_alignment": str(output.resolve()), "analysis_index": str(final_index.resolve()),
            "normalization_reference": reference_fingerprint(reference), "source_inputs": inputs,
            "alignment_policy": policy, "bwa_mem2_index": bwa_mem2_index,
            "library_provenance": library_provenance,
            "duplicate_marking": duplicate_marking,
            "source_file_count": len(rows),
            "source_runs": [{"run_id": row["run_id"], "format": row["format"], "path": row["path"],
                             "read_group_id": row.get("read_group_id")} for row in rows],
            "source_record_tag_counts": source_stats, "output_record_tag_counts": output_stats,
            "tools": tools, "source_io_execution": source_io_execution,
            "threads": {"align": align_threads, "sort": sort_threads, "extract": extract_threads},
            "metrics": metrics, "log": str(log_path),
        }
        write_json_atomic(qc_json, payload)
        if keep_run_bams:
            keep = output.parent / f"{sample}.stage04_runs"
            if keep.exists():
                shutil.rmtree(keep)
            shutil.move(str(temp_dir), keep)
        else:
            shutil.rmtree(temp_dir)
        return payload
    except BaseException as exc:
        failure = {"schema": "mfm_stage04_alignment_qc_v4", "sample_id": sample,
                   "source_sample_id": source_sample_id, "status": "failed",
                   "canonical_format": "bam", "error": str(exc), "temporary_directory": str(temp_dir),
                   "alignment_policy": policy, "duplicate_marking": duplicate_plan,
                   "library_provenance": library_provenance, "log": str(log_path)}
        write_json_atomic(qc_json, failure)
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-json", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--qc-json", required=True, type=Path)
    parser.add_argument("--minimap2", required=True)
    parser.add_argument("--bwa-mem2", default="bwa-mem2")
    parser.add_argument("--illumina-aligner", choices=("bwa-mem2", "minimap2"), default="bwa-mem2")
    parser.add_argument("--bwa-mem2-index-root", type=Path)
    parser.add_argument("--picard", default="picard")
    parser.add_argument("--illumina-duplicate-policy", choices=("auto", "mark", "preserve"), default="auto")
    parser.add_argument("--markdup-heap-gb", type=int, default=16)
    parser.add_argument("--samtools", required=True)
    parser.add_argument("--align-threads", type=int, default=32)
    parser.add_argument("--sort-threads", type=int, default=8)
    parser.add_argument("--sort-memory", default="2G")
    parser.add_argument("--extract-threads", type=int, default=4)
    parser.add_argument("--min-quality-values", type=int, default=2)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--keep-run-bams", action="store_true")
    args = parser.parse_args()
    contract = json.loads(args.contract_json.read_text())
    payload = prepare(contract, output=args.output, qc_json=args.qc_json,
                      minimap2=args.minimap2, bwa_mem2=args.bwa_mem2,
                      illumina_aligner=args.illumina_aligner,
                      bwa_mem2_index_root=args.bwa_mem2_index_root, picard=args.picard,
                      illumina_duplicate_policy=args.illumina_duplicate_policy,
                      markdup_heap_gb=args.markdup_heap_gb, samtools=args.samtools,
                      align_threads=args.align_threads, sort_threads=args.sort_threads,
                      sort_memory=args.sort_memory, extract_threads=args.extract_threads,
                      min_quality_values=args.min_quality_values,
                      overwrite=args.overwrite, keep_run_bams=args.keep_run_bams)
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
