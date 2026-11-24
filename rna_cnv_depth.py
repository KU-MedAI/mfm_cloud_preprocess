"""HAP×염색체 bedcov 병렬 실행. 통계의 모집단·필터·합산 순서는 바꾸지 않는다."""
from __future__ import annotations

import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent/"src"))

from rna_resources import worker_count


def _bedcov(task: dict, samtools: str, stop, processes, lock) -> list[tuple[str, float]]:
    start_time = time.monotonic()
    label = f"HAP{task['hap'] + 1}:{task['contig']}:{task.get('part',0)}"
    print(f"[cnv:start] {label} intervals={len(task['rows']):,}", flush=True)
    with tempfile.TemporaryDirectory(prefix="rna_cnv_") as temp:
        bed = Path(temp) / "genes.bed"
        with bed.open("w") as handle:
            for contig, start, end, gene in task["rows"]:
                handle.write(f"{contig}\t{start}\t{end}\t{gene}\n")
        command = [samtools, "bedcov"]
        if task["reference"]:
            command += ["--reference", task["reference"]]
        command += [str(bed), str(task["bam"])]
        if task.get('native_tool'):
            command=[task['native_tool'],str(task['bam']),str(task['reference'] or '-'),task['contig'],str(bed),str(task.get('decode_threads',0))]
        # native 출력은 파일로 받고 종료 후 파싱해 계산과 Python PIPE 복사의 CPU 중첩을 막음.
        with (Path(temp)/'counts.tsv').open('w+') as output:
            with lock:
                if stop.is_set():
                    raise RuntimeError(f"{label}: 다른 task 실패로 취소")
                proc = subprocess.Popen(command, stdout=output if task.get("native_tool") else subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                processes.add(proc)
            try:
                stdout, stderr = proc.communicate()
                if task.get('native_tool'):output.seek(0);stdout=output.read()
            finally:
                with lock:
                    processes.discard(proc)
        if proc.returncode:
            raise RuntimeError(f"{label}: samtools bedcov exit {proc.returncode}\n{stderr[-2000:]}")
        covered = []
        for line in stdout.splitlines():
            fields = line.split("\t")
            if len(fields) >= 5:
                covered.append((fields[3], float(fields[4])))
        # samtools가 일부 행만 반환해도 완료로 처리하지 않는다.
        if len(covered) != len(task["rows"]):
            raise RuntimeError(f"{label}: bedcov 행 수 불일치 {len(covered)} != {len(task['rows'])}")
    print(f"[cnv:done] {label} elapsed={time.monotonic() - start_time:.1f}s", flush=True)
    return covered


def _native_bedcov(task,samtools):
    return _bedcov(task,samtools,threading.Event(),set(),threading.Lock())


def _native_task(task,samtools):
    from integrate_wgs_rna.rna_task_chunks import cpu_seconds,timing
    started=time.monotonic();cpu=cpu_seconds()
    rows=_native_bedcov(task,samtools)
    return rows,timing(started,cpu,hap=task['hap'],contig=task['contig'],part=task['part'],rows=len(rows))


def collect_depth(inputs: list[dict], samtools: str, requested: str, *, execution=None) -> list[dict[str, float]]:
    tasks = []
    for hap, data in enumerate(inputs):
        by_contig: dict[str, list] = {}
        for row in data["rows"]:
            by_contig.setdefault(row[0], []).append(row)
        for contig, rows in by_contig.items():
            tasks.append({"hap": hap, "contig": contig, "rows": rows,
                          "bam": data["bam"], "reference": data["reference"]})
    if not tasks:
        raise ValueError("bedcov 입력 구간이 없습니다")
    import os
    from integrate_wgs_rna.rna_wgs_native import select,plan
    from integrate_wgs_rna.cpu_executor import cpu_results
    tool=select('RNA_CNV_BACKEND','rna_cnv_sum')
    results={}
    if tool:
        check=subprocess.run([samtools,'quickcheck',*[str(data['bam']) for data in inputs]],capture_output=True,text=True)
        if check.returncode:raise RuntimeError('CNV 입력 quickcheck 실패: '+check.stderr[-2000:])
        from integrate_wgs_rna.rna_task_chunks import partition_rows,task_layout
        work=[]
        for original,task in enumerate(tasks):
            for part,indices in enumerate(partition_rows(task['rows'],lambda row:(row[1],row[2]))):
                work.append(dict(task,rows=[task['rows'][i] for i in indices],indices=indices,original=original,part=part))
        layout=task_layout();budget=plan(len(work),requested=requested,decode_threads=0 if layout=='chunks' else None)
        for task in work:task.update(native_tool=tool,decode_threads=budget['decode_threads'])
        ordered=[[None]*len(task['rows']) for task in tasks];traces=[]
        for index,value,error in cpu_results(_native_task,[(task,samtools) for task in work],budget['workers_initial'],elastic=True,
                max_workers=budget['workers_upper'],cpus_per_worker=budget['cpu_per_worker'],memory_per_worker_bytes=budget['worker_address_space_bytes'],parent_cpus=budget['parent_cpus']):
            if error:raise error
            values,trace=value;task=work[index];traces.append(trace)
            print(f'[cnv:progress] {len(traces)}/{len(work)} completed [native]',flush=True)
            for position,row in zip(task['indices'],values):ordered[task['original']][position]=row
        if any(row is None for rows in ordered for row in rows):raise ValueError('CNV 분할 결과 행 누락')
        results=dict(enumerate(ordered))
        if execution is not None:execution.update(layout=layout,budget=budget,tasks=traces)
        import json
        print('[RNA-CNV-TASKS] '+json.dumps(dict(layout=layout,budget=budget,tasks=traces)),flush=True)
    else:
        jobs = worker_count(requested, len(tasks))
        print(f"[cnv] tasks={len(tasks)} workers={jobs}; HAP 전체 중앙값 정책 유지", flush=True)
        results = {}
        stop, lock, processes = threading.Event(), threading.Lock(), set()
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(_bedcov, task, samtools, stop, processes, lock): i
                       for i, task in enumerate(tasks)}
            try:
                for done, future in enumerate(as_completed(futures), 1):
                    results[futures[future]] = future.result()
                    print(f"[cnv:progress] {done}/{len(tasks)} completed", flush=True)
            except BaseException:
                # 비정상 task가 나오면 이 함수가 생성한 다른 bedcov도 종료한다.
                with lock:
                    stop.set()
                    for proc in processes:
                        if proc.poll() is None:
                            proc.terminate()
                for future in futures:
                    future.cancel()
                raise
    sums: list[dict[str, float]] = [{} for _ in inputs]
    # 완료 순서가 아닌 원래 BED 정렬 순서로 float 덧셈. 여러 chr에 걸친 gene도 동일하게 집계.
    for i, task in enumerate(tasks):
        dest = sums[task["hap"]]
        for gene, value in results[i]:
            dest[gene] = dest.get(gene, 0.0) + value
    return sums
