#!/usr/bin/env python3
"""17 — Picard/Parabricks 중복 표시 (duplicate 표시만, 삭제하지 않음) + 최종 HAP BAM.

    python 17_mark_duplicates.py --sample S --hap 1 --align-json <W>/rna/hap1.align.json \
        --out-dir <R>/HAP1 --work-dir <W>/rna --log-dir <W>/logs

입력   16 의 STAR BAM(hap{h}.align.json)
출력   <align_dir>/Aligned.sortedByCoord.out.md.bam(.bai), <align_dir>/<S>.HAP{h}.markduplicates.metrics.txt
       <out-dir>/<S>.rna.personal.HAP{h}.bam(.bai)   최종 BAM (21 입력)
       <work-dir>/hap{h}.markdup.json                 metrics 요약 (19 가 rna_realign.json 에 넣는다)

옵션은 GTEx `rnaseq/src/run_MarkDuplicates.py` 39–48줄과 같다: PROGRAM_RECORD_ID=null MAX_RECORDS_IN_RAM=500000
SORTING_COLLECTION_SIZE_RATIO=0.25 ASSUME_SORT_ORDER=coordinate TAGGING_POLICY=DontTag OPTICAL_DUPLICATE_PIXEL_DISTANCE=100.
duplicate 는 flag/metrics 만 남기고 지우지 않는다(GTEx 도 REMOVE_DUPLICATES 미지정 = false).
duplicate 정책: dense coverage(18)는 GTEx bam2coverage.py 처럼 duplicate 를 **포함**하고,
ASE count(21)는 ASE QC 목적으로 duplicate를 **제외**한다(한 HAP에서라도 dup인 fragment를 양쪽에서 제외).
★ -Djava.io.tmpdir 는 TMP_DIR 와 다르다. TMP_DIR 은 정렬 spill, java.io.tmpdir 은 JVM 이 snappy 같은 native .so 를
  푸는 곳이다. 이 호스트의 /tmp 는 noexec 라 후자를 옮기지 않으면 UnsatisfiedLinkError 로 죽는다.
"""
from __future__ import annotations

import argparse
import json
import time
import os
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from rna_cache import cache_fields, cache_valid, file_fingerprint, read_cache, tool_fingerprint

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rna_common import die, log, nonempty, read_json, run, write_json  # noqa: E402
import rna_profile  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from integrate_wgs_rna.markdup_resources import markdup_plan
from integrate_wgs_rna import parabricks_bam, gpu_rna


def metric_values(path: Path) -> dict[str, str]:
    """Picard metrics 의 첫 표(헤더 + 첫 행)를 dict 로."""
    header = None
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if header is None:
            header = f
            continue
        return dict(zip(header, f))
    return {}


def num(v: str | None) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def link_or_copy_atomic(source: Path, destination: Path) -> bool:
    """Keep both legacy paths; avoid duplicate storage when hard links are available."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    part = Path(str(destination) + ".part")
    part.unlink(missing_ok=True)
    try:
        os.link(source, part)
        linked = True
    except OSError:
        shutil.copyfile(source, part)
        linked = False
    os.replace(part, destination)
    return linked


def primary_record_counts(path: Path, threads: int = 8) -> tuple[int, int]:
    """Count primary and primary-duplicate records in one BAM pass.

    This is equivalent to ``samtools view -c -F 2308`` and the same query with
    ``-f 1024``: unmapped, secondary, and supplementary records are excluded.
    QC-fail records remain included, matching the original command.
    """
    from integrate_wgs_rna.record_stats_parallel import native_tool,NativeUnavailable
    import subprocess
    try:tool=native_tool('rna_primary')
    except NativeUnavailable:tool=None
    if tool:
        total=int(os.environ.get('MFM_SAMPLE_CPU_BUDGET',str(threads+2)))
        result=subprocess.run([str(tool),str(path),str(min(max(0,total-2),max(0,threads)))],capture_output=True,text=True)
        if result.returncode:raise RuntimeError('RNA 독립 flag QC 실패: '+result.stderr[-2000:])
        return tuple(json.loads(result.stdout))
    primary = duplicates = 0
    with pysam.AlignmentFile(str(path), "rb", threads=max(1, threads)) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            primary += 1
            duplicates += int(read.is_duplicate)
    return duplicates, primary


def index_and_counts(path, samtools, execution):
    """안정된 BAM의 색인과 독립 QC를 같은 CPU grant 안에서 겹침."""
    cpu = execution['cpu_grant']
    if cpu < 5:
        auxiliary=max(0,cpu-2)
        run([samtools,'index','-@',str(auxiliary),str(path)],what='MarkDuplicates BAM index')
        return primary_record_counts(path,threads=auxiliary), dict(parallel=False,cpu_grant=cpu)
    children=cpu-1
    index_cpus=children//2;qc_cpus=children-index_cpus
    plan=dict(parallel=True,cpu_grant=cpu,parent_cpus=1,index_cpus=index_cpus,qc_cpus=qc_cpus)
    log('MARKDUP-INDEX-QC-PLAN '+json.dumps(plan,sort_keys=True))
    with ThreadPoolExecutor(max_workers=2) as pool:
        index=pool.submit(run,[samtools,'index','-@',str(max(0,index_cpus-2)),str(path)],
                          what='MarkDuplicates BAM index')
        counts=pool.submit(primary_record_counts,path,threads=max(0,qc_cpus-2))
        index.result()
        return counts.result(),plan


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--hap", required=True, choices=("1", "2"))
    ap.add_argument("--align-json", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--log-dir", required=True, type=Path)
    ap.add_argument("--picard", default=os.environ.get("PICARD_BIN", "picard"))
    ap.add_argument("--heap-gb", type=int, default=int(os.environ.get("RNA_MARKDUP_HEAP_GB", "8")))
    ap.add_argument("--samtools", default=os.environ.get("SAMTOOLS", "samtools"))
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    S, h, st = a.sample, a.hap, a.samtools
    use_parabricks = parabricks_bam.enabled()
    if not use_parabricks and not (os.path.isfile(a.picard) and os.access(a.picard, os.X_OK)) and shutil.which(a.picard) is None:
        die(f"{S} HAP{h} Picard 없음: {a.picard} (conda env update -n jj_mfm -f environment.yml 또는 PICARD_BIN 지정 필요)")
    if a.heap_gb < 1:
        die(f"RNA_MARKDUP_HEAP_GB는 양의 정수여야 함: {a.heap_gb}")
    al = read_json(a.align_json)
    ad, star_bam = Path(al["align_dir"]), Path(al["star_bam"])
    if not nonempty(star_bam):
        die(f"{S} HAP{h} STAR BAM 없음: {star_bam} (16 먼저)")
    if use_parabricks:parabricks_bam.check_markdup_validation()
    integrated=al.get('integrated_markdup') if use_parabricks else None
    integrated_metrics=gpu_rna.validate_integrated_markdup(integrated,star_bam) if integrated is not None else None
    execution = markdup_plan(a.heap_gb,star_bam)
    started=time.monotonic()
    log("MARKDUP-PLAN " + json.dumps(execution,sort_keys=True))
    md = ad / "Aligned.sortedByCoord.out.md.bam"
    mt = ad / f"{S}.HAP{h}.markduplicates.metrics.txt"
    inputs = {"stage": 17, "sample": S, "hap": h,
              "star_bam": file_fingerprint(star_bam),
              "markdup_args": list(rna_profile.MARKDUP_ARGS),
              "bam_backend": "parabricks" if use_parabricks else "cpu",
              "integrated_markdup": integrated,
              "backend_script": file_fingerprint(parabricks_bam.__file__, content=True),
              "picard": tool_fingerprint(a.picard), "samtools": tool_fingerprint(st),
              "script": file_fingerprint(__file__, content=True)}
    inputs['bam_compression_level']=None if use_parabricks else 1
    record_path = a.work_dir / f"hap{h}.markdup.json"
    regenerate = a.force or not nonempty(md, Path(str(md) + ".bai"), mt) or not cache_valid(read_cache(record_path), inputs)
    backend_execution = read_cache(record_path).get("backend_execution") if not regenerate else None
    verified_counts=None;index_qc_plan=None
    if regenerate:
        log(f"  HAP{h}: {'Parabricks markdup' if use_parabricks else 'Picard MarkDuplicates'}")
        tmpd = a.work_dir / f"picard_tmp_hap{h}"; tmpd.mkdir(parents=True, exist_ok=True)
        part_bam, part_mt = Path(str(md) + ".part.bam"), Path(str(mt) + ".part.txt")
        if integrated_metrics is not None:
            log(f"  HAP{h}: STAR 통합 duplicate 결과 사용; 복사·index·독립 flag QC 수행")
            # source inode를 link하면 align cache의 ctime이 바뀌므로 독립 복사함.
            shutil.copyfile(star_bam,part_bam)
            shutil.copyfile(integrated_metrics,part_mt)
            backend_execution={'backend':'parabricks','integrated_with_alignment':True,
                               'alignment_json':str(a.align_json),'remove_duplicates':False,
                               'baseline_sort_order':'coordinate','optical_duplicate_pixel_distance':100}
        elif use_parabricks:
            reference=parabricks_bam.alignment_reference(al, tmpd)
            backend_execution=parabricks_bam.mark_duplicates(star_bam, reference, part_bam, part_mt, tmpd, a.log_dir / f"parabricks.markdup.HAP{h}.log")
        else:
            run([a.picard, f"-Xmx{a.heap_gb}g", f"-XX:ActiveProcessorCount={execution['java_active_processors']}", f"-Djava.io.tmpdir={tmpd}", "MarkDuplicates",
                 f"I={star_bam}", f"O={part_bam}", f"TMP_DIR={tmpd}", f"M={part_mt}",
                 *rna_profile.MARKDUP_ARGS,"COMPRESSION_LEVEL=1"],
                log_path=a.log_dir / f"17.markduplicates.HAP{h}.log", what=f"{S} HAP{h} Picard MarkDuplicates")
        verified_counts,index_qc_plan=index_and_counts(part_bam,st,execution)
        os.replace(part_bam, md); os.replace(str(part_bam) + ".bai", str(md) + ".bai"); os.replace(part_mt, mt)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    final = a.out_dir / f"{S}.rna.personal.HAP{h}.bam"
    link_or_copy_atomic(md, final)
    link_or_copy_atomic(Path(str(md) + ".bai"), Path(str(final) + ".bai"))

    m = metric_values(mt)
    dup_primary, primary = (verified_counts if verified_counts is not None
                            else primary_record_counts(md,threads=max(0,execution['cpu_grant']-2)))
    pd, pu = num(m.get("READ_PAIR_DUPLICATES")), num(m.get("UNPAIRED_READ_DUPLICATES"))
    write_json(record_path, {
        "tool": "Parabricks markdup" if use_parabricks else "Picard MarkDuplicates",
        "bam_backend": "parabricks" if use_parabricks else "cpu",
        "backend_execution": backend_execution,
        "resource_execution": {**execution,"wall_seconds":time.monotonic()-started,
                               "cache_reused":not regenerate,"index_qc_plan":index_qc_plan},
        "remove_duplicates": False,
        "marked_bam": str(md),
        "final_bam": str(final),
        "metrics_path": str(mt),
        "read_pairs_examined": int(num(m.get("READ_PAIRS_EXAMINED"))),
        "read_pair_duplicates": int(pd),
        "unpaired_reads_examined": int(num(m.get("UNPAIRED_READS_EXAMINED"))),
        "unpaired_read_duplicates": int(pu),
        "duplicate_read_records_from_metrics": int(2 * pd + pu),
        "primary_duplicate_records": dup_primary,
        "primary_records": primary,
        "percent_duplication": num(m.get("PERCENT_DUPLICATION")),
        **cache_fields(inputs, [md, Path(str(md) + ".bai"), mt]),
    })
    log(f"  HAP{h}: duplicate {num(m.get('PERCENT_DUPLICATION')):.4f} → {final}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
