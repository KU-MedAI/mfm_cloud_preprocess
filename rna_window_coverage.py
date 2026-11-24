#!/usr/bin/env python3
"""전장 RNA coverage에서 WGS retained-window coverage를 파생한다.

새 depth를 계산하지 않는다. 같은 HAP의 coverage.zarr를 source로 읽고
WGS windows.tsv의 kept 행 순서대로 lossless subset을 만든다.
"""
from __future__ import annotations

import hashlib
import os
import shutil
import uuid
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc,blosc

import rna_profile

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from integrate_wgs_rna.window_quality import depth_metrics
from integrate_wgs_rna.window_contract import read_windows  # noqa: E402


def sha256_file(path: Path, chunk: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def retained_windows(path: Path, hap: str, lengths: dict[str, int],
                     contigs: list[str]) -> list[dict]:
    allowed = set(contigs)
    rows: list[dict] = []
    for raw in read_windows(path):
        if raw["haplotype"] != hap or raw["contig"] not in allowed:
            continue
        if raw["kept"].strip().lower() not in {"1", "true"}:
            continue
        start, end = int(raw["personal_start"]), int(raw["personal_end"])
        if raw["contig"] not in lengths or end > lengths[raw["contig"]]:
            raise ValueError(f"RNA window가 HAP FASTA 범위를 벗어남: {raw}")
        rows.append({
            "contig": raw["contig"],
            "haplotype": hap,
            "window_index": int(raw["window_index"]),
            "personal_start": start,
            "personal_end": end,
            "ref_start": int(raw.get("ref_start") or -1),
            "ref_end": int(raw.get("ref_end") or -1),
        })
    if not rows:
        raise ValueError(f"{hap} retained WGS window가 없음: {path}")
    keys = [
        (row["contig"], row["window_index"], row["personal_start"], row["personal_end"])
        for row in rows
    ]
    if len(keys) != len(set(keys)):
        raise ValueError(f"{hap} retained WGS window key 중복")
    observed = {row["contig"] for row in rows}
    if observed != allowed:
        raise ValueError(
            f"{hap} retained WGS contig 불일치 expected={sorted(allowed)} observed={sorted(observed)}"
        )
    return rows


def _text_array(group: zarr.Group, name: str, values: list[str], width: int) -> None:
    data = np.asarray(values, dtype=f"<U{width}")
    group.create_dataset(
        name, data=data, chunks=(min(max(1, len(data)), 4096),), overwrite=True
    )


def _commit_directory(part: Path, output: Path) -> None:
    """완성된 임시 directory를 기존 출력과 같은 filesystem에서 교체한다."""
    previous = output.with_name(f".{output.name}.previous.{os.getpid()}.{uuid.uuid4().hex}")
    moved_previous = False
    try:
        if output.exists() or output.is_symlink():
            if output.is_symlink() or not output.is_dir():
                raise ValueError(f"window coverage output이 정상 directory가 아님: {output}")
            os.replace(output, previous)
            moved_previous = True
        os.replace(part, output)
    except BaseException:
        if moved_previous and previous.exists() and not output.exists():
            os.replace(previous, output)
        raise
    finally:
        if previous.exists():
            shutil.rmtree(previous)


def write_window_coverage(
    *,
    source_path: Path,
    windows_path: Path,
    output_path: Path,
    sample: str,
    hap: str,
    lengths: dict[str, int],
    contigs: list[str],
    chunk: int,
    power: float,
    clip: float,
    source_buffers: Path | None = None,
) -> dict:
    """전장 source를 window별로 slice하여 flat depth + 좌표 계약을 기록한다."""
    if hap not in {"HAP1", "HAP2"}:
        raise ValueError(f"잘못된 HAP label: {hap}")
    source_path = source_path.absolute()
    windows_path = windows_path.absolute()
    output_path = output_path.absolute()
    if source_path == output_path or source_path in output_path.parents or output_path in source_path.parents:
        raise ValueError("전장 coverage와 window coverage 경로는 독립적이어야 함")

    source = zarr.open_group(str(source_path), mode="r")
    rna_profile.validate_coverage_metadata(source.attrs)
    if str(source.attrs.get("sample")) != sample or str(source.attrs.get("hap")) != hap:
        raise ValueError("전장 coverage sample/HAP 불일치")
    if source.attrs.get("coordinate_space") != "individual (personal haplotype)":
        raise ValueError("전장 coverage 좌표계 불일치")
    if int(source.attrs.get("total_depth", 0)) <= 0:
        raise ValueError("전장 coverage가 완료되지 않았거나 비어 있음")
    for contig in contigs:
        if contig not in source:
            raise ValueError(f"전장 coverage contig 누락: {contig}")
        if source[contig].shape != (lengths[contig],) or source[contig].dtype != np.dtype(np.uint32):
            raise ValueError(f"전장 coverage shape/dtype 불일치: {contig}")

    rows = retained_windows(windows_path, hap, lengths, contigs)
    offsets = [0]
    for row in rows:
        offsets.append(offsets[-1] + row["personal_end"] - row["personal_start"])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    part = output_path.with_name(f".{output_path.name}.part.{os.getpid()}.{uuid.uuid4().hex}")
    if part.exists():
        raise RuntimeError(f"임시 window coverage 경로 충돌: {part}")

    from integrate_wgs_rna.cpu_executor import WorkerBudget
    from integrate_wgs_rna.gpu_compare_resources import cpu_limit
    available=cpu_limit(int(os.environ.get('MFM_SAMPLE_CPU_BUDGET',os.environ.get('MFM_CPU_BUDGET','68'))))
    # contig worker 종료 후 반환된 CPU를 압축 전에 다시 빌림. 작은 chunk는 추가 thread 불필요함.
    requested=min(available,max(1,chunk*4//(1<<20))+1)
    budget=WorkerBudget(requested,requested,True)
    if source_buffers:budget.base_memory+=sum(path.stat().st_size for path in Path(source_buffers).glob('*.bin'))
    grant=0
    while grant<1:
        grant=budget.capacity(requested)
        if grant<1:__import__('time').sleep(.2)
    previous_threads=blosc.get_nthreads();blosc.set_nthreads(max(1,grant-1))
    print(f'[cov:window-write] cpu_upper={grant} blosc_threads={max(1,grant-1)} chunk_bytes={chunk*4}',flush=True)
    compressor = Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE)
    try:
        root = zarr.open_group(str(part), mode="w")
        total_bases = offsets[-1]
        depth = root.zeros(
            "depth", shape=(total_bases,), dtype=np.uint32,
            chunks=(min(chunk, max(1, total_bases)),), compressor=compressor, overwrite=True,
        )
        root.create_dataset(
            "offsets", data=np.asarray(offsets, dtype=np.int64),
            chunks=(min(len(offsets), 4096),), compressor=compressor, overwrite=True,
        )
        for name in ("window_index", "personal_start", "personal_end", "ref_start", "ref_end"):
            root.create_dataset(
                name, data=np.asarray([row[name] for row in rows], dtype=np.int64),
                chunks=(min(len(rows), 4096),), compressor=compressor, overwrite=True,
            )
        _text_array(root, "contig", [row["contig"] for row in rows], 32)
        _text_array(root, "haplotype", [hap] * len(rows), 8)
        _text_array(
            root,
            "window_key",
            [
                f'{row["contig"]}:{hap}:{row["window_index"]}:'
                f'{row["personal_start"]}-{row["personal_end"]}'
                for row in rows
            ],
            96,
        )

        # flat array chunk를 한 번만 압축함. 1Mb window마다4Mb chunk를 재압축하지 않음.
        pending=np.empty(min(chunk,max(1,total_bases)),dtype=np.uint32)
        pending_n=written=0
        total_depth = nonzero_positions = nonzero_sum = 0
        quality = []
        depth_digest = hashlib.sha256()
        buffers={};validated=set();current_contig=None;current_buffer=None
        if source_buffers:
            for entry in Path(source_buffers).glob('*.json'):
                record=__import__('json').loads(entry.read_text());buffers[record['contig']]=record
        for index, row in enumerate(rows):
            contig=row['contig']
            if contig in buffers:
                record=buffers[contig];path=Path(record['path']);stamp=path.stat()
                if stamp.st_size!=lengths[contig]*4 or stamp.st_size!=record['bytes'] or stamp.st_mtime_ns!=record['mtime_ns']:raise ValueError('RNA buffer 변경/크기 오류')
                if contig not in validated:
                    if sha256_file(path)!=source.attrs['depth_sha256_by_contig'][contig]:raise ValueError('RNA buffer SHA 불일치')
                    validated.add(contig)
                if current_contig!=contig:
                    current_buffer=np.memmap(path,mode='r',dtype='<u4',shape=(lengths[contig],));current_contig=contig
                values=current_buffer[row['personal_start']:row['personal_end']]
            else:
                values = source[contig][row['personal_start']:row['personal_end']]
            begin, end = offsets[index], offsets[index + 1]
            pos=0
            while pos<len(values):
                size=min(len(values)-pos,len(pending)-pending_n)
                pending[pending_n:pending_n+size]=values[pos:pos+size]
                pos+=size;pending_n+=size
                if pending_n==len(pending):
                    depth[written:written+pending_n]=pending
                    written+=pending_n;pending_n=0
            quality.append(depth_metrics(values, end-begin, backend="cpu"))
            nonzero = values > 0
            total_depth += int(values.sum())
            nonzero_positions += int(nonzero.sum())
            nonzero_sum += int(values[nonzero].sum())
            depth_digest.update(memoryview(np.asarray(values, dtype="<u4")).cast("B"))
        if pending_n:depth[written:written+pending_n]=pending[:pending_n]
        if total_depth <= 0 or nonzero_positions <= 0:
            raise ValueError(f"{sample} {hap}: retained WGS windows의 RNA coverage가 0")

        for name in ("depth_mean", "depth_zero_fraction", "depth_max"):
            root.create_dataset(name, data=np.asarray([row[name] for row in quality],
                                dtype=np.uint32 if name == "depth_max" else np.float64),
                                compressor=compressor, chunks=(min(len(rows),4096),), overwrite=True)
        root.attrs.update({
            "schema": "mfm_rna_window_coverage_v1",
            "sample": sample,
            "hap": hap,
            "source_coverage_zarr": str(source_path),
            "source_bam_signature": source.attrs.get("source_bam_signature"),
            "source_fai_signature": source.attrs.get("source_fai_signature"),
            "coordinate_space": "individual_personal_haplotype",
            "window_contract": "exact Stage WGS windows.tsv row order filtered by haplotype and kept=true",
            "windows_tsv": str(windows_path),
            "windows_tsv_sha256": sha256_file(windows_path),
            "window_count": len(rows),
            "retained_bases": total_bases,
            "total_depth": total_depth,
            "nonzero_positions": nonzero_positions,
            "nonzero_mean": nonzero_sum / nonzero_positions,
            "nonzero_mean_scope": "retained_wgs_windows",
            "depth_sha256": depth_digest.hexdigest(),
            "depth_sha256_encoding": "flat retained-window raw little-endian uint32 C-order",
            "norm_power": power,
            "norm_clip_at": clip,
            **rna_profile.coverage_metadata(),
        })
        _commit_directory(part, output_path)
    except BaseException:
        if part.exists():
            shutil.rmtree(part)
        raise
    finally:
        blosc.set_nthreads(previous_threads)
        budget.base_memory=4*(1<<30)
        budget.capacity(0,grow=False)

    return {
        "window_count": len(rows),
        "retained_bases": offsets[-1],
        "total_depth": total_depth,
        "nonzero_positions": nonzero_positions,
    }
