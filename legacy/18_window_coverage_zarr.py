#!/usr/bin/env python3
"""18 — WGS windows.tsv와 1:1인 HAP별 RNA coverage Zarr를 생성한다.

전장 chromosome 배열을 최종 보존하지 않는다. WGS Stage 11에서 선택한 개인 HAP
window만 읽으며, windows.tsv의 HAP별 행 순서가 곧 RNA depth 행 순서다.
"""
from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rna_profile  # noqa: E402
from rna_resources import file_signature  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from integrate_wgs_rna.window_contract import read_windows  # noqa: E402


def sha256(path: Path, chunk: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def load_windows(path: Path, hap: str, lengths: dict[str, int], contigs: list[str]) -> list[dict]:
    allowed = set(contigs)
    rows: list[dict] = []
    for raw in read_windows(path):
        if raw["haplotype"] != hap or raw["contig"] not in allowed:
            continue
        # dropped window는 audit에만 남긴다. 모델 track은 compact와 같은 kept 행만 사용한다.
        if str(raw["kept"]).strip().lower() not in {"1", "true"}:
            continue
        start, end = int(raw["personal_start"]), int(raw["personal_end"])
        if end > lengths.get(raw["contig"], 0):
            raise ValueError(f"invalid window: {raw}")
        rows.append({
            "contig": raw["contig"], "haplotype": hap,
            "window_index": int(raw["window_index"]),
            "personal_start": start, "personal_end": end,
            "ref_start": int(raw.get("ref_start") or -1),
            "ref_end": int(raw.get("ref_end") or -1),
        })
    if not rows:
        raise ValueError(f"no retained WGS windows for {hap}: {path}")
    keys = [(r["contig"], r["window_index"], r["personal_start"], r["personal_end"]) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError(f"duplicate WGS window keys for {hap}")
    return rows


def write_text_array(group: zarr.Group, name: str, values: list[str], width: int) -> None:
    data = np.asarray(values, dtype=f"<U{width}")
    group.create_dataset(name, data=data, chunks=(min(max(1, len(data)), 4096),), overwrite=True)


def coverage_source(path: Path, args, lengths: dict, contigs: list) -> zarr.Group:
    """같은 BAM·HAP·필터로 계산한 전장 raw depth만 재사용한다."""
    src, dst = path.resolve(), args.out.resolve()
    if src == dst or src in dst.parents or dst in src.parents:
        raise ValueError("coverage source와 output은 겹치지 않는 독립 경로여야 합니다")
    source = zarr.open_group(str(path), mode="r")
    rna_profile.validate_coverage_metadata(source.attrs)
    expected = {
        "sample": args.sample, "hap": args.hap,
        "coordinate_space": "individual (personal haplotype)",
        "min_mapq_applied": 255,
        "mate_overlap_policy": "counted twice (samtools depth -s not used)",
        "deletion_policy": "counted as covered (samtools depth -J)",
        "source_bam_signature": file_signature(args.bam),
        "source_fai_signature": file_signature(args.fai),
        "samtools_depth_args": list(rna_profile.COVERAGE_DEPTH_ARGS),
    }
    for key, value in expected.items():
        if source.attrs.get(key) != value:
            raise ValueError(f"coverage source 불일치: {key}")
    if source.attrs.get("total_depth", 0) <= 0:
        raise ValueError("coverage source가 완료되지 않았거나 비어 있습니다")
    for contig in contigs:
        if (contig not in source or source[contig].shape != (lengths[contig],)
                or source[contig].dtype != np.dtype(np.uint32)):
            raise ValueError(f"coverage source contig shape/dtype 불일치: {contig}")
    return source


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--fai", required=True, type=Path)
    ap.add_argument("--windows-tsv", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--chunk", type=int, default=4194304)
    ap.add_argument("--power", type=float, default=0.75)
    ap.add_argument("--clip", type=float, default=10.0)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--hap", required=True, choices=["HAP1", "HAP2"])
    ap.add_argument("--samtools", default=os.environ.get("SAMTOOLS", "samtools"))
    ap.add_argument("--coverage-source-zarr", type=Path,
                    help="동일 BAM에서 만든 전장 coverage. 생략하면 독립 samtools depth 계산")
    a = ap.parse_args()

    lengths: dict[str, int] = {}
    with a.fai.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            lengths[fields[0]] = int(fields[1])
    contigs = [c for c in a.contigs.replace(",", " ").split() if c]
    missing_contigs = sorted(set(contigs) - set(lengths))
    if missing_contigs:
        raise SystemExit(f"[FAIL] FASTA FAI missing contigs: {','.join(missing_contigs)}")

    rows = load_windows(a.windows_tsv, a.hap, lengths, contigs)
    # 입력 계약 검사는 출력 store를 열어 덮어쓰기 전에 끝낸다.
    source_signature = file_signature(a.bam)
    source = coverage_source(a.coverage_source_zarr, a, lengths, contigs) if a.coverage_source_zarr else None
    by_contig: dict[str, list[tuple[int, dict]]] = defaultdict(list)
    offsets = [0]
    for index, row in enumerate(rows):
        by_contig[row["contig"]].append((index, row))
        offsets.append(offsets[-1] + row["personal_end"] - row["personal_start"])

    comp = Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE)
    root = zarr.open(str(a.out), mode="w")
    total_bases = offsets[-1]
    depth = root.zeros(
        "depth", shape=(total_bases,), dtype=np.uint32,
        chunks=(min(a.chunk, max(1, total_bases)),), compressor=comp, overwrite=True,
    )
    root.create_dataset("offsets", data=np.asarray(offsets, dtype=np.int64),
                        chunks=(min(len(offsets), 4096),), compressor=comp, overwrite=True)
    for name in ("window_index", "personal_start", "personal_end", "ref_start", "ref_end"):
        root.create_dataset(
            name, data=np.asarray([row[name] for row in rows], dtype=np.int64),
            chunks=(min(len(rows), 4096),), compressor=comp, overwrite=True,
        )
    write_text_array(root, "contig", [row["contig"] for row in rows], 32)
    write_text_array(root, "haplotype", [a.hap] * len(rows), 8)
    write_text_array(
        root, "window_key",
        [f'{row["contig"]}:{a.hap}:{row["window_index"]}:{row["personal_start"]}-{row["personal_end"]}'
         for row in rows], 96,
    )

    depth_args = [arg for arg in rna_profile.COVERAGE_DEPTH_ARGS if arg != "-a"]
    nz_sum = nz_n = total_depth = 0
    for contig in contigs:
        indexed = by_contig.get(contig, [])
        if not indexed:
            continue
        indexed.sort(key=lambda item: item[1]["personal_start"])
        starts = [row["personal_start"] for _, row in indexed]
        ends = [row["personal_end"] for _, row in indexed]
        local_offsets = [0]
        for start, end in zip(starts, ends):
            local_offsets.append(local_offsets[-1] + end - start)
        local = np.zeros(local_offsets[-1], dtype=np.uint32)

        if source is not None:
            # RAM에 contig 한 개를 읽어 작은 window별 디스크 재독해를 피한다.
            chromosome = source[contig][:]
            for i, (start, end) in enumerate(zip(starts, ends)):
                local[local_offsets[i]:local_offsets[i + 1]] = chromosome[start:end]
            del chromosome
        else:
            a.out.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=f".{a.hap}.{contig}.bed", dir=a.out.parent, delete=False
            ) as bed:
                bed_path = Path(bed.name)
                for start, end in zip(starts, ends):
                    bed.write(f"{contig}\t{start}\t{end}\n")
            try:
                command = [
                    a.samtools, "depth", *depth_args, "-b", str(bed_path), "-r", contig, str(a.bam)
                ]
                with tempfile.TemporaryFile(mode="w+t") as stderr:
                    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=stderr, text=True)
                    try:
                        assert proc.stdout is not None
                        for line in proc.stdout:
                            fields = line.rstrip("\n").split("\t")
                            pos0, value = int(fields[1]) - 1, int(fields[2])
                            i = bisect.bisect_right(starts, pos0) - 1
                            if i >= 0 and pos0 < ends[i]:
                                local[local_offsets[i] + pos0 - starts[i]] = value
                        rc = proc.wait()
                        if rc:
                            stderr.seek(0)
                            raise RuntimeError(f"samtools depth failed ({rc}) {contig}: {stderr.read()[-2000:]}")
                    finally:
                        if proc.poll() is None:
                            proc.terminate()
                            proc.wait()
                        if proc.stdout is not None:
                            proc.stdout.close()
            finally:
                bed_path.unlink(missing_ok=True)

        for local_i, (global_i, _) in enumerate(indexed):
            begin, end = local_offsets[local_i], local_offsets[local_i + 1]
            depth[offsets[global_i]:offsets[global_i + 1]] = local[begin:end]
        mask = local > 0
        nz_sum += int(local[mask].sum())
        nz_n += int(mask.sum())
        total_depth += int(local.sum())
        print(f"  {contig}: windows {len(indexed):,} / depth {int(local.sum()):,} / >0 {int(mask.sum()):,}")

    if total_depth <= 0 or nz_n <= 0:
        raise SystemExit(f"[FAIL] {a.sample} {a.hap}: retained WGS windows의 RNA coverage가 0입니다")
    if file_signature(a.bam) != source_signature:
        raise RuntimeError("window coverage 계산 중 source BAM이 변경됐습니다")
    nonzero_mean = nz_sum / nz_n
    root.attrs.update({
        "schema": "mfm_rna_window_coverage_v1",
        "sample": a.sample, "hap": a.hap,
        "source_bam": str(a.bam.resolve()),
        "coordinate_space": "individual_personal_haplotype",
        "window_contract": "exact row order after filtering Stage 11 windows.tsv by haplotype and kept=true",
        "windows_tsv": str(a.windows_tsv.resolve()),
        "windows_tsv_sha256": sha256(a.windows_tsv),
        "window_count": len(rows), "retained_bases": total_bases,
        "total_depth": total_depth, "nonzero_positions": nz_n,
        "nonzero_mean": nonzero_mean, "nonzero_mean_scope": "retained_wgs_windows",
        "norm_power": a.power, "norm_clip_at": a.clip,
        **rna_profile.coverage_metadata(),
    })
    print(json.dumps({
        "sample": a.sample, "hap": a.hap, "windows": len(rows),
        "retained_bases": total_bases, "total_depth": total_depth,
        "nonzero_positions": nz_n, "nonzero_mean": nonzero_mean,
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
