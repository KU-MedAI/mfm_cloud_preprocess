#!/usr/bin/env python3
"""개인 유전체 영구 자산(HAP FASTA·chain·QC) manifest 작성과 검증.

    python tools/personal_assets.py manifest --individual-fasta-dir <final>/<S>/individual_fasta \
        --contigs "chr1 … chr22" [--reference-fasta <HAP 을 만든 reference>]
    python tools/personal_assets.py verify --individual-fasta-dir <…>

★ 왜 영구 보존하는가 (2026-09-28)
  09 는 SNV+indel 을 적용하므로 HAP 좌표는 reference 와 다르다. chain 이 있어야
  - GENCODE 가 바뀌어도 Stage 14에서 GTF만 다시 lift 할 수 있고,
  - WGS 01–12 를 다시 돌리지 않고 RNA 13 이후만 재실행할 수 있다.
  그래서 WGS Stage 24/RNA Stage 26이 호출하는 공용 cleanup backend는 이 파일들을
  절대 지우지 않고, 지우기 전에 이 manifest로 온전함을 확인한다.

영구 자산(PERMANENT)
  HAP1|HAP2/<chr>.fa, <chr>.fa.fai, <chr>.chain, personal_fasta_qc.json, permanent_assets.json(이 manifest)
manifest 에는 파일별 크기·sha256, chain 헤더(tName/tSize/qName/qSize), chain 길이 일치 검사,
reference fingerprint(personal_fasta_qc.json 값, 없으면 여기서 계산)를 적는다.
"""
from __future__ import annotations

import argparse
import os
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

MANIFEST = "permanent_assets.json"
SCHEMA = "mfm_personal_permanent_assets_v1"


def sha256(path: Path, chunk: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def read_fai(path: Path) -> dict[str, int]:
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in path.open() if l.strip()}


def chain_header(path: Path) -> dict:
    f = path.open().readline().split()
    if len(f) < 12 or f[0] != "chain":
        raise ValueError(f"chain 헤더가 아니다: {path}")
    return {"tName": f[2], "tSize": int(f[3]), "qName": f[7], "qSize": int(f[8])}


def reference_fingerprint(reference: Path) -> dict:
    """src/integrate_wgs_rna/personal_fasta.py reference_fingerprint 와 같은 정의(FASTA stat + FAI sha256)."""
    reference = reference.resolve()
    fai = Path(str(reference) + ".fai")
    st = reference.stat()
    return {"path": str(reference), "size_bytes": st.st_size, "mtime_ns": st.st_mtime_ns,
            "fai_path": str(fai.resolve()), "fai_size_bytes": fai.stat().st_size, "fai_sha256": sha256(fai)}


def file_record(path: Path) -> dict:
    return {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def permanent_paths(root: Path, contigs: list[str]) -> list[Path]:
    out = [root / "personal_fasta_qc.json"]
    for hap in ("HAP1", "HAP2"):
        for c in contigs:
            out += [root / hap / f"{c}.fa", root / hap / f"{c}.fa.fai", root / hap / f"{c}.chain"]
    return out


def build(root: Path, contigs: list[str], reference: Path | None, *, workers: int = 1) -> dict:
    errors: list[str] = []
    qc_path = root / "personal_fasta_qc.json"
    qc = json.loads(qc_path.read_text()) if qc_path.is_file() else {}
    if not qc:
        errors.append(f"personal_fasta_qc.json 없음: {qc_path}")
    ref = reference or (Path(qc["reference_fasta"]) if qc.get("reference_fasta") else None)
    ref_fp = qc.get("reference_fingerprint")
    ref_len: dict[str, int] = {}
    if ref is not None and Path(str(ref) + ".fai").is_file():
        ref_len = read_fai(Path(str(ref) + ".fai"))
        if not ref_fp:
            ref_fp = reference_fingerprint(Path(ref))
    else:
        errors.append(f"reference .fai 를 찾지 못함: {ref}")
    files, chains = {}, {}
    paths=permanent_paths(root,contigs)
    valid=[]
    for p in paths:
        if not p.is_file() or p.stat().st_size==0:
            errors.append(f"영구 자산 없음 또는 비어 있음: {p}");continue
        valid.append(p)
    with ThreadPoolExecutor(max_workers=max(1,min(workers,len(os.sched_getaffinity(0))))) as pool:
        for p,record in zip(valid,pool.map(file_record,valid)):
            files[str(p.relative_to(root))]=record
    for hap in ("HAP1", "HAP2"):
        for c in contigs:
            ch, fai = root / hap / f"{c}.chain", root / hap / f"{c}.fa.fai"
            if not ch.is_file() or not fai.is_file():
                continue
            h = chain_header(ch)
            hap_len = read_fai(fai).get(c)
            check = {**h, "tName_ok": h["tName"] == c, "qName_ok": h["qName"] == c,
                     "qSize_eq_hap_fai": h["qSize"] == hap_len,
                     "tSize_eq_reference": (h["tSize"] == ref_len.get(c)) if ref_len else None}
            for k in ("tName_ok", "qName_ok", "qSize_eq_hap_fai", "tSize_eq_reference"):
                if check[k] is False:
                    errors.append(f"{hap} {c} chain 검사 실패: {k} ({h})")
            chains[f"{hap}/{c}"] = check
    return {"schema": SCHEMA, "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "individual_fasta_dir": str(root), "contigs": contigs, "sequence_variant_policy": qc.get("sequence_variant_policy"),
            "reference_fasta": str(ref) if ref else None, "reference_fingerprint": ref_fp,
            "policy": "PERMANENT: cleanup backend는 이 파일들을 지우지 않는다. RNA Stage 13 이후 재실행과 GTF 재-lift에 쓴다.",
            "files": files, "chains": chains, "status": "pass" if not errors else "fail", "errors": errors}


def verify(root: Path, *, workers: int = 1) -> dict:
    man = json.loads((root / MANIFEST).read_text())
    bad = []
    def check(item):
        rel,rec=item;p=root/rel
        if not p.is_file():return f"없음: {rel}"
        if p.stat().st_size!=rec['size_bytes'] or sha256(p)!=rec['sha256']:return f"변경됨: {rel}"
        return None
    with ThreadPoolExecutor(max_workers=max(1,min(workers,len(os.sched_getaffinity(0))))) as pool:
        bad=[error for error in pool.map(check,man['files'].items()) if error is not None]
    return {"status": "pass" if not bad and man.get("status") == "pass" else "fail", "problems": bad,
            "files_checked": len(man["files"])}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("manifest", "verify"))
    ap.add_argument("--individual-fasta-dir", required=True, type=Path)
    ap.add_argument("--contigs", default="")
    ap.add_argument("--reference-fasta", type=Path)
    a = ap.parse_args()
    root = a.individual_fasta_dir
    if a.command == "manifest":
        contigs = a.contigs.split()
        if not contigs:
            print("--contigs 필요", file=sys.stderr); return 2
        man = build(root, contigs, a.reference_fasta)
        tmp = root / (MANIFEST + ".part")
        tmp.write_text(json.dumps(man, ensure_ascii=False, indent=2) + "\n")
        tmp.replace(root / MANIFEST)
        print(f"[personal-assets] {root}: {len(man['files'])} files, status={man['status']}"
              + ("" if not man["errors"] else f" errors={man['errors'][:3]}"))
        return 0 if man["status"] == "pass" else 1
    res = verify(root)
    print(json.dumps(res, ensure_ascii=False))
    return 0 if res["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
