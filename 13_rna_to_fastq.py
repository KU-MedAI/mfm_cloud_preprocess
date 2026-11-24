#!/usr/bin/env python3
"""13 — Stage 02 입력 계약을 확인하고 RNA BAM/CRAM에서 paired FASTQ를 복원한다.

    python 13_rna_to_fastq.py --sample S --rna <RNA BAM/CRAM> --out-dir <W>/rna --log-dir <W>/logs \
        [--decode-reference <CRAM 원 reference>] [--sm-expected "S,public_id"] [--force]

입력   RNA BAM/CRAM (manifest rna_path; 이미 어딘가에 정렬된 RNA — GEUVADIS/GTEx 는 GRCh38 STAR BAM)
출력   <out-dir>/R1.fq, R2.fq (비압축), <out-dir>/fastq_qc.json (schema mfm_rna_fastq_qc_v2)
재사용 R1·R2·fastq_qc.json 이 모두 있으면 다시 만들지 않는다(입력 검사는 매번 한다). AS2-14 도 이 FASTQ 를 재사용한다.

⓪ 입력 계약 재확인
  - samtools quickcheck, @RG SM, paired fraction 전수 검사는 Stage 02 RNA preflight가 담당한다.
  - Stage 13은 preflight JSON의 경로·크기·mtime·inode를 확인해 검사 뒤 바뀐 입력을 거부한다.
  - CRAM 디코딩에는 원 정렬 reference와 .fai가 필요하다.
① BAM → FASTQ (paired 보존)
  - 원본은 coordinate-sorted 다. `samtools collate` 로 이름 순서를 맞추지 않으면 페어가 깨진다.
  - `samtools fastq -n` : R1/R2 는 짝이 있는 read 만, singleton·unpaired 는 세고 버린다.
    secondary/supplementary 는 samtools fastq 기본 동작으로 빠진다(read 당 1회).
  - 과거 파일럿 메모에는 14분 대 2.3분이 기록되어 있으나 입력·서버·측정일 정보가 없어
    현재 성능 기준으로 사용하지 않는다.
  - paired 비율은 기본적으로 QC에만 기록한다. strict 호환 모드에서만 과거 cutoff를 적용한다.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from rna_cache import cache_fields, cache_valid, file_fingerprint, optional_file, read_cache, tool_fingerprint

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rna_common import die, log, nonempty, write_json  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from integrate_wgs_rna import parabricks_bam


def fastq_backend():
    requested=os.environ.get('RNA_FASTQ_BACKEND','auto')
    if requested not in {'auto','samtools','parabricks'}:raise ValueError('RNA_FASTQ_BACKEND는 auto|samtools|parabricks임')
    return 'samtools' if requested=='auto' else requested


def samtools_cpu_plan(collate,fastq):
    limits=parabricks_bam.resource_plan()
    total=limits['cpus']
    if total<2:raise ValueError('paired FASTQ pipe에는 최소2CPU가 필요함')
    extra=max(0,total-4)
    a=min(max(0,collate),(extra+1)//2);b=min(max(0,fastq),extra-a)
    peak=int((4*1024**3+(a+b)*64*1024**2)*1.1)
    if peak>limits['host_memory_bytes']:raise ValueError('FASTQ pipe 예상RAM이예약예산보다큼')
    return {'cpu_grant':total,'collate_aux_threads':a,'fastq_aux_threads':b,
            'cpu_upper':2+a+b+int(a>0)+int(b>0),
            'predicted_peak_ram_bytes':peak,'ram_budget_bytes':limits['host_memory_bytes']}


def fastq_stats(path: Path) -> tuple[int, int]:
    """(레코드 수, 최대 read 길이). 4 줄 단위가 아니면 형식 오류로 중단."""
    if not nonempty(path):
        return 0, 0
    out = subprocess.run(["awk", "NR%4==2 && length($0)>m{m=length($0)} END{if (NR%4) exit 2; print NR/4, m+0}",
                          str(path)], capture_output=True, text=True)
    if out.returncode != 0:
        die(f"FASTQ 레코드 형식 오류: {path}")
    n, m = out.stdout.split()
    return int(n), int(m)


def prepare_missing_preflight(a):
    """WGS read QC를생략한경로에서도기존RNA입력검증전체를수행함."""
    if os.environ.get('MFM_WGS_PAYLOAD')!='sequence_only' or os.environ.get('MFM_RNA_OUTPUT_MODE')!='coverage_only':return False
    if not a.gtf or not a.contigs:raise ValueError('RNA 전용 preflight에 GTF/contigs가 필요함')
    cpus=parabricks_bam.resource_plan()['cpus']
    command=[sys.executable,str(Path(__file__).with_name('rna_preflight.py')),
             '--sample',a.sample,'--expected-samples',a.sm_expected or a.sample,'--rna',str(a.rna),
             '--gtf',str(a.gtf),'--contigs',a.contigs,'--samtools',a.samtools,
             '--hap-read-mode',a.hap_read_mode,'--threads',str(max(1,cpus-1)),
             '--min-paired-fraction',str(a.min_paired_fraction),'--ratio-qc-policy',a.ratio_qc_policy,
             '--sm-placeholders',a.sm_placeholders,'--skip-hap-assets',
             '--require-sm-match' if a.require_sm_match else '--no-require-sm-match',
             '--output',str(a.preflight_json)]
    if a.decode_reference:command+=['--decode-reference',a.decode_reference]
    subprocess.run(command,check=True)
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--rna", required=True)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--log-dir", required=True, type=Path)
    ap.add_argument("--decode-reference", default="")
    ap.add_argument("--sm-expected", default="", help="쉼표 구분. 비우면 --sample")
    ap.add_argument("--sm-placeholders", default=os.environ.get("RNA_SM_PLACEHOLDERS", "sm1"))
    ap.add_argument("--require-sm-match", type=int, default=int(os.environ.get("RNA_REQUIRE_SM_MATCH", "1")))
    ap.add_argument("--min-paired-fraction", type=float, default=float(os.environ.get("RNA_MIN_PAIRED_FRACTION", "0.90")))
    ap.add_argument(
        "--ratio-qc-policy", choices=("record_only", "strict"),
        default=os.environ.get("QC_RATIO_POLICY", "record_only"),
    )
    ap.add_argument("--samtools", default=os.environ.get("SAMTOOLS", "samtools"))
    ap.add_argument("--gtf", type=Path)
    ap.add_argument("--contigs", default="")
    ap.add_argument("--hap-read-mode", default="all_to_both")
    ap.add_argument("--preflight-json", required=True,
                    help="Stage 02 rna_preflight.json; source가 바뀌면 FASTQ 복원을 거부한다")
    ap.add_argument("--collate-threads", type=int, default=8)
    ap.add_argument("--fastq-threads", type=int, default=4)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    S, rna, st = a.sample, a.rna, a.samtools
    if not rna or rna == "-" or not Path(rna).is_file() or Path(rna).stat().st_size == 0:
        die(f"{S} RNA BAM 없음: {rna}")
    od, ld = a.out_dir, a.log_dir
    od.mkdir(parents=True, exist_ok=True); ld.mkdir(parents=True, exist_ok=True)

    # ── ⓪ Stage 02 입력 계약 재검증 ──
    # quickcheck·header·paired fraction 전수 검사는 Stage 02에서 한 번 수행한다.
    # 여기서는 같은 inode/size/mtime의 source인지 확인해 바뀐 입력의 캐시 재사용을 막는다.
    contract_path = Path(a.preflight_json)
    if not nonempty(contract_path) and not prepare_missing_preflight(a):
        die(f"{S} Stage 02 RNA preflight 없음: {contract_path}")
    contract = json.loads(contract_path.read_text())
    source = contract.get("rna", {})
    stat = Path(rna).resolve().stat()
    if (contract.get("status") != "pass"
            or Path(source.get("path", "")).resolve() != Path(rna).resolve()
            or int(source.get("size", -1)) != stat.st_size
            or int(source.get("mtime_ns", -1)) != stat.st_mtime_ns
            or int(source.get("inode", -1)) != stat.st_ino):
        die(f"{S} Stage 02 이후 RNA 입력이 변경됨: {rna}")
    ref_args: list[str] = []
    if rna.endswith(".cram"):
        dr = a.decode_reference
        if not (dr and nonempty(Path(dr)) and nonempty(Path(dr + ".fai"))):
            die(f"{S} RNA 입력이 CRAM 인데 decode reference(.fai 포함)가 없음: "
                "manifest rna_reference_fasta 또는 RNA_DECODE_REFERENCE 지정 필요")
        ref_args = ["--reference", dr]
    sm = {
        "status": contract.get("checks", {}).get("sm_status", "unknown"),
        "observed": contract.get("checks", {}).get("read_group_samples", []),
        "expected_any_of": contract.get("checks", {}).get("expected_samples", []),
        "source": str(contract_path),
    }

    # ── ① BAM → FASTQ ──
    r1, r2, qc = od / "R1.fq", od / "R2.fq", od / "fastq_qc.json"
    decode_reference = a.decode_reference if ref_args else ""
    selected_backend=fastq_backend()
    inputs = {"stage": 13, "sample": S, "source": file_fingerprint(rna),
              "preflight": {"schema": contract.get("schema"), "rna": source,
                            "checks": contract.get("checks", {})},
              "decode_reference": optional_file(decode_reference),
              "decode_fai": optional_file(decode_reference + ".fai" if decode_reference else None, content=True),
              "min_paired_fraction": a.min_paired_fraction, "ratio_qc_policy": a.ratio_qc_policy,
              "bam_backend": "parabricks" if parabricks_bam.enabled() else "cpu",
              "fastq_backend": selected_backend, "fastq_backend_requested": os.environ.get("RNA_FASTQ_BACKEND","auto"),
              "backend_script": file_fingerprint(parabricks_bam.__file__, content=True),
              "samtools": tool_fingerprint(st), "script": file_fingerprint(__file__, content=True)}
    if not a.force and nonempty(r1, r2, qc) and cache_valid(read_cache(qc), inputs):
        log(f"{S}: FASTQ 재사용 ({od})")
        return 0
    log(f"{S}: RNA BAM → FASTQ")
    parts = {k: od / f"{n}.part" for k, n in
             (("r1", "R1.fq"), ("r2", "R2.fq"), ("un", "discard_unpaired.fq"), ("sg", "discard_singleton.fq"))}
    for p in parts.values():
        p.unlink(missing_ok=True)
    backend_execution = None
    if selected_backend=='parabricks':
        backend_execution = parabricks_bam.rna_fastq(rna, decode_reference, parts, ld / "parabricks.bam2fq.log")
    else:
        cpu=samtools_cpu_plan(a.collate_threads,a.fastq_threads)
        from integrate_wgs_rna.rna_fastq_native import reserve_conversion
        cpu['conversion_reservation']=reserve_conversion(cpu)
        log('FASTQ-PLAN '+json.dumps(cpu,sort_keys=True));started=time.monotonic()
        with (ld / "13.collate.err").open("w") as e1, (ld / "13.fastq.err").open("w") as e2:
            col = subprocess.Popen([st, "collate", "-u", "-O", "-@", str(cpu["collate_aux_threads"]), *ref_args, rna,
                                    str(od / "collate_tmp")], stdout=subprocess.PIPE, stderr=e1)
            fq = subprocess.Popen([st, "fastq", "-@", str(cpu["fastq_aux_threads"]), "-1", str(parts["r1"]), "-2", str(parts["r2"]),
                                   "-0", str(parts["un"]), "-s", str(parts["sg"]), "-n", "-"],
                                  stdin=col.stdout, stderr=e2)
            col.stdout.close()
            rc_fq, rc_col = fq.wait(), col.wait()
        if rc_col != 0 or rc_fq != 0:
            die(f"{S} FASTQ 추출 실패 (collate {rc_col}, fastq {rc_fq})")
        backend_execution={'backend':'samtools','compression':'uncompressed','resource_plan':cpu,
                           'elapsed_seconds':time.monotonic()-started,'gpu_kernels_requested':False}
    stats=backend_execution.get('fastq_stats') if backend_execution else None
    if stats is None:
        from integrate_wgs_rna.rna_fastq_native import batch
        cpus=parabricks_bam.resource_plan()['cpus']
        stats=dict(zip(['r1','r2','sg','un'],batch([(parts[k],None,None) for k in ['r1','r2','sg','un']],cpus)))
    p1,m1=stats['r1']['records'],stats['r1']['max_read_length'];p2,m2=stats['r2']['records'],stats['r2']['max_read_length']
    sg,un=stats['sg']['records'],stats['un']['records']
    if p1 != p2 or p1 <= 0:
        die(f"{S} FASTQ pair 계약 위반: R1={p1} R2={p2}")
    tot = 2 * p1 + sg + un
    paired_fraction = round(2 * p1 / tot if tot else 0.0, 6)
    below_reference = paired_fraction < a.min_paired_fraction
    if a.ratio_qc_policy == "strict" and below_reference:
        die(f"{S} RNA paired 비율 부족: {paired_fraction} < {a.min_paired_fraction} (single-end RNA는 지원하지 않음)")
    os.replace(parts["r1"], r1); os.replace(parts["r2"], r2)
    write_json(qc, {
        "schema": "mfm_rna_fastq_qc_v2",
        "input_alignment": rna,
        "bam_backend": "parabricks" if parabricks_bam.enabled() else "cpu",
              "fastq_backend": selected_backend, "fastq_backend_requested": os.environ.get("RNA_FASTQ_BACKEND","auto"),
        "backend_execution": backend_execution,
        "rna_preflight_contract": str(contract_path),
        "rna_preflight_contract_sha256": contract.get("contract_sha256"),
        "input_decode_reference": a.decode_reference or None,
        "input_sm_check": sm,
        "r1_fastq": str(r1),
        "r2_fastq": str(r2),
        "r1_pairs": p1,
        "r2_pairs": p2,
        "singleton_discard_records": sg,
        "unpaired_discard_records": un,
        "paired_fraction": paired_fraction,
        "paired_fraction_reference": a.min_paired_fraction,
        "paired_fraction_below_reference": below_reference,
        "paired_fraction_qc_action": a.ratio_qc_policy,
        "max_read_length": max(m1, m2),
        "pair_contract_pass": p1 == p2,
        **cache_fields(inputs, [r1, r2]),
    })
    parts["un"].unlink(missing_ok=True); parts["sg"].unlink(missing_ok=True)
    for p in od.glob("collate_tmp*.bam"):
        p.unlink(missing_ok=True)
    log(f"{S}: FASTQ {p1:,} pair, paired {paired_fraction:.4f}, max read {max(m1, m2)} bp")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
