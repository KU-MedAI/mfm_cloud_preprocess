"""염색체 단위 WhatsHap 계산을 병렬화하고 기존 순서로 결합한다.

염색체 전체를 처리하여 PS/linked-read/supplementary 계산 범위를 보존한다.
기존 --skip-missing-contigs 동작을 따르며 unplaced read는 한번만 복사한다.
실패한 shard와 기존 부분 출력은 보존하고 성공한 결과만 게시한다.
"""
from __future__ import annotations

import gzip
import fcntl
import json
import os
import shutil
import subprocess
import time
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pysam

from .haplotag_resources import reserve_tail


TAG_HEADER = b"#readname\thaplotype\tphaseset\tchromosome\n"


class ChromosomeBudget:
    """공유 CPU 풀의 동일 임대를 잠금 아래 증액·반환한다."""

    def __init__(self, cpus: int, count: int):
        value = os.environ.get("MFM_CPU_POOL_LEASE_FILE")
        self.path = Path(value) if value else None
        self.remaining = count
        self.initial = cpus
        self.active = self.peak = 0
        self.lock = threading.Lock()

    def capacity(self, grow: bool) -> int:
        if self.path is None:
            return min(self.initial, max(1, self.remaining))
        root = self.path.parent.parent
        with (root / "pool.lock").open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            row = json.loads(self.path.read_text())
            capacity = min(int(row["cpus"]), max(1, self.remaining))
            if row.get('cpu_backfill'):
                capacity = min(capacity, self.initial, int(row['cpu_backfill']['cpus']))
            # 대기 중인 다른 단계가 있으면 그 단계가 먼저 자원을 받는다.
            if grow and not row.get('cpu_backfill') and not any((root / "waiters").glob("*.json")):
                other = 0
                for path in (root / "leases").glob("*.json"):
                    if path != self.path:
                        try:
                            other += int(json.loads(path.read_text())["cpus"])
                        except FileNotFoundError:
                            pass
                free_for_us = max(0, int(row["total_cpus"]) - other)
                capacity = max(capacity, min(self.remaining, free_for_us))
            row["cpus"] = capacity
            row["peak_cpus"] = max(int(row.get("peak_cpus", 0)), capacity)
            candidate = self.path.with_suffix("." + uuid.uuid4().hex + ".tmp")
            candidate.write_text(json.dumps(row) + "\n")
            os.replace(candidate, self.path)
            return capacity

    def acquire(self):
        while True:
            with self.lock:
                if self.active < self.capacity(grow=True):
                    self.active += 1
                    self.peak = max(self.peak, self.active)
                    return
            time.sleep(0.1)

    def release(self):
        with self.lock:
            self.active -= 1
            self.remaining -= 1
            self.capacity(grow=False)


def run_parallel(config, sample: str, alignment: Path, tag_vcf: Path,
                 reference: Path, output_bam: Path, taglist: Path,
                 contract: dict, cpus: int) -> dict:
    """각 subprocess는 계산/압축을 한 코어에서 수행하여 합산 예산을 지킨다."""
    cpus = max(1, int(cpus))
    work = output_bam.parent / "parallel_haplotag" / (
        time.strftime("%Y%m%dT%H%M%S") + "_" + uuid.uuid4().hex[:12])
    work.mkdir(parents=True, exist_ok=False)
    report_path = work / "execution.json"
    with pysam.VariantFile(str(tag_vcf)) as variants:
        vcf_contigs = set(variants.header.contigs)
    with pysam.AlignmentFile(str(alignment), reference_filename=str(reference)) as source:
        statistics = source.get_index_statistics()
        populated = {s.contig for s in statistics if s.total > 0}
        # CRAM index는 contig별 통계를 제공하지 않을 수 있어 전체 header를 사용한다.
        if alignment.suffix.lower() == ".cram":
            populated = set(source.references)
        contigs = [c for c in source.references if c in populated and c in vcf_contigs]
        source_header = source.header.to_dict()
    jobs = min(cpus, max(1, len(contigs)))
    budget = ChromosomeBudget(cpus, len(contigs))
    report = {"schema": "mfm_parallel_haplotag_v1", "sample": sample,
              "input_alignment": str(alignment), "cpus": cpus, "jobs": jobs,
              "contigs": contigs, "status": "running", "started": time.time(),
              "tag_supplementary": config.tag_supplementary,
              "ignore_read_groups": contract["ignore_read_groups"]}

    def save() -> None:
        report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")

    save()

    def shard(item):
        budget.acquire()
        try:
            return run_shard(item)
        finally:
            budget.release()

    def run_shard(item):
        i, chrom = item
        bam = work / f"{i:04d}.bam"
        tags = work / f"{i:04d}.tags.tsv.gz"
        cmd = [config.whatshap_bin, "haplotag", "--sample", sample,
               "--reference", str(reference), "--output-threads", "0",
               "--skip-missing-contigs", "--regions", chrom,
               "--output-haplotag-list", str(tags), "-o", str(bam)]
        if config.tag_supplementary:
            cmd.append("--tag-supplementary=copy-primary")
        if contract["ignore_read_groups"]:
            cmd.append("--ignore-read-groups")
        cmd += [str(tag_vcf), str(alignment)]
        env = os.environ.copy()
        # 외부 gzip worker 없이 각 shard를 한 CPU에서 처리한다.
        # WhatsHap은 gzip 이외 확장자에서는 일반 텍스트를 출력한다.
        plain_tags = tags.with_suffix("")
        cmd[cmd.index("--output-haplotag-list") + 1] = str(plain_tags)
        with (work / f"{i:04d}.log").open("wb") as log:
            subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([config.samtools_bin, "quickcheck", str(bam)], check=True)
        return bam, plain_tags

    try:
        # 공유 임대가 있는 경우 여유 코어가 생기면 미실행 염색체를 즉시 추가한다.
        max_workers = max(1, len(contigs)) if budget.path else jobs
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            shards = list(pool.map(shard, enumerate(contigs)))
        reserve_tail(config, sample, "combine_and_tags")
        # WhatsHap의 fetch('*')와 같은 unplaced read를 태그 변경 없이 한번 복사한다.
        unmapped = work / "unmapped.bam"
        # 모든 염색체가 끝났으면 임대는 1CPU로 축소되어 다른 계산에 반환된다.
        subprocess.run([config.samtools_bin, "view", "-b", "-@", "0",
                        "-T", str(reference), "-o", str(unmapped), str(alignment), "*"], check=True)
        bams = [bam for bam, _ in shards] + [unmapped]
        header = work / "header.sam"
        with pysam.AlignmentFile(str(bams[0])) as first:
            assert first.header.to_dict().get("SQ") == source_header.get("SQ")
            assert first.header.to_dict().get("RG") == source_header.get("RG")
            header.write_text(str(first.header))
        candidate = work / "combined.bam"
        subprocess.run([config.samtools_bin, "cat", "-h", str(header),
                        "-o", str(candidate), *map(str, bams)], check=True)
        subprocess.run([config.samtools_bin, "quickcheck", str(candidate)], check=True)
        candidate_tags = work / "combined.tags.tsv.gz"
        with gzip.open(candidate_tags, "wb", compresslevel=1) as target:
            target.write(TAG_HEADER)
            for _, tags in shards:
                with tags.open("rb") as src:
                    if src.readline() != TAG_HEADER:
                        raise ValueError(f"태그 목록 header 불일치: {tags}")
                    shutil.copyfileobj(src, target, 16 << 20)
        # 이전 부분 태그 목록은 이름/경로를 바꾸지 않고 복제해 보존한다.
        if taglist.exists():
            os.link(taglist, work / "previous.tags.tsv.gz")
        # 같은 파일시스템의 완성 후보만 최종 이름으로 게시한다.
        os.replace(candidate_tags, taglist)
        os.replace(candidate, output_bam)
        report.update(status="success", finished=time.time(), peak_concurrent=budget.peak,
                      output_bam=str(output_bam), taglist=str(taglist))
        save()
        return report
    except BaseException as exc:
        report.update(status="failed", finished=time.time(), error=str(exc))
        save()
        raise
