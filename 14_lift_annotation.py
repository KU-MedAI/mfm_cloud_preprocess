#!/usr/bin/env python3
"""14 — GENCODE GTF 를 chain 으로 HAP 개인 좌표로 옮긴다.

    python 14_lift_annotation.py --chain-dir <individual_fasta/HAP1> --gtf <gencode.v47.gtf.gz> \
        --out-gtf <HAP1/annotation.gtf> --out-exons <HAP1/exons.bygene.tsv> --qc-json <HAP1/gtf_lift_qc.json> \
        [--partner-chain-dir <individual_fasta/HAP2>] [--reference-fai <ref.fa.fai>] [--contigs "chr1 … chr22"] [--reuse]

왜: STAR index의 annotated splice junction(15)과 allele count의 exon→gene 배정(21)이 HAP 좌표여야 한다.
    reference 좌표 GTF 를 그대로 쓰면 HAP 의 indel 만큼 어긋난다. RNA read 자체는 GTF 와 무관하게 HAP FASTA 전체에 정렬된다.

chain: `bcftools consensus --chain`(09) 출력. target = reference, query = HAP.
결실에 걸린 feature: 정렬되는 염기만 남긴다(시작 = 첫 정렬 염기, 끝 = 마지막 정렬 염기, `lift_trimmed "1"`).
    결실된 염기는 그 HAP 에 없으므로 잘린 exon 이 그 HAP 의 실제 exon 이다. feature 전체가 결실이면 그 줄만 뺀다.
좌표계 검사(어긋나면 중단): GTF 헤더가 GRCh38 이외 build 를 명시, chain tName·qName ≠ contig,
    chain tSize ≠ --reference-fai 길이, chain qSize ≠ HAP .fai 길이, 두 HAP chain 의 tSize 불일치, feature end > tSize.
재사용(--reuse): QC JSON 이 같은 schema·입력 GTF 이고 출력이 있으면 다시 만들지 않는다.
"""
from __future__ import annotations
import argparse, bisect, gzip, json, os, re, sys
from pathlib import Path
from rna_cache import cache_fields, cache_valid, file_fingerprint, optional_file, read_cache

SCHEMA = "mfm_gtf_lift_qc_v5"


def sha256_file(path: str, chunk: int = 8 << 20) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()
NON_GRCH38_BUILD = re.compile(r"GRCh37|hg19|hg18|CHM13|T2T|hs1\b", re.IGNORECASE)


def load_chain(path: Path):
    """((t_starts, q_starts, sizes), header) — reference→HAP 정렬 블록과 chain 헤더."""
    ts, qs, sz = [], [], []
    with open(path) as fh:
        hdr = fh.readline().split()
        if len(hdr) < 12 or hdr[0] != "chain":
            raise ValueError(f"{path}: chain 헤더가 아니다")
        header = {"tName": hdr[2], "tSize": int(hdr[3]), "qName": hdr[7], "qSize": int(hdr[8])}
        t, q = int(hdr[5]), int(hdr[10])
        for line in fh:
            p = line.split()
            if not p:
                break
            if len(p) == 3:
                s, dt, dq = int(p[0]), int(p[1]), int(p[2])
                ts.append(t); qs.append(q); sz.append(s)
                t += s + dt; q += s + dq
            else:
                s = int(p[0]); ts.append(t); qs.append(q); sz.append(s); break
    return (ts, qs, sz), header


def map_interval(s0: int, e0: int, ts, qs, sz):
    """reference 0-based 닫힌 구간 [s0, e0] → (HAP 시작, HAP 끝(닫힘), 정렬 염기 수, 경계가 잘렸는지). 정렬 염기가 없으면 None."""
    i = max(bisect.bisect_right(ts, s0) - 1, 0)
    if ts[i] + sz[i] <= s0:
        i += 1
    first = last = ref_first = ref_last = None
    aligned = 0
    while i < len(ts) and ts[i] <= e0:
        lo, hi = max(s0, ts[i]), min(e0, ts[i] + sz[i] - 1)
        if lo <= hi:
            aligned += hi - lo + 1
            if first is None:
                first, ref_first = qs[i] + (lo - ts[i]), lo
            last, ref_last = qs[i] + (hi - ts[i]), hi
        i += 1
    if first is None:
        return None
    return first, last, aligned, (ref_first != s0 or ref_last != e0)


def read_fai(path: Path) -> dict[str, int]:
    return {f[0]: int(f[1]) for f in (l.split("\t") for l in open(path)) if len(f) >= 2}


def load_chain_set(chain_dir: Path, contigs: list[str], label: str, ref_lengths: dict[str, int] | None):
    """contig 별 chain 을 읽고 좌표계를 검사한다. (maps, checks)."""
    maps, checks, missing = {}, [], []
    for c in contigs:
        p = chain_dir / f"{c}.chain"
        if not p.exists():
            missing.append(c); continue
        m, h = load_chain(p)
        if h["tName"] != c or h["qName"] != c:
            raise ValueError(f"{label} {p}: chain 이름이 contig 와 다르다 (tName={h['tName']} qName={h['qName']})")
        if ref_lengths is not None:
            if c not in ref_lengths:
                raise ValueError(f"--reference-fai 에 {c} 가 없다")
            if ref_lengths[c] != h["tSize"]:
                raise ValueError(f"{label} {c}: chain tSize {h['tSize']} != reference 길이 {ref_lengths[c]}"
                                 " — HAP FASTA 를 만든 reference 와 다르다")
        fai = chain_dir / f"{c}.fa.fai"
        q_checked = fai.is_file()
        if q_checked and read_fai(fai).get(c) != h["qSize"]:
            raise ValueError(f"{label} {c}: chain qSize {h['qSize']} != HAP FASTA 길이 {read_fai(fai).get(c)}")
        maps[c] = (m, h)
        checks.append({"hap": label, "contig": c, "tSize": h["tSize"], "qSize": h["qSize"],
                       "tSize_vs_reference": "pass" if ref_lengths is not None else "not_checked",
                       "qSize_vs_hap_fai": "pass" if q_checked else "not_checked"})
    if missing:
        raise FileNotFoundError(f"GTF lift에 필요한 {label} chain이 없습니다: " + ",".join(missing))
    return maps, checks


def get_attr(attrs: str, key: str) -> str:
    for kv in attrs.split(";"):
        kv = kv.strip()
        if kv.startswith(key + " "):
            parts = kv.split('"')
            return parts[1] if len(parts) > 1 else kv.split()[1]
    return ""


def gtf_build_check(gtf: str, opener) -> dict:
    header = []
    with opener(gtf, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            header.append(line.rstrip("\n"))
    text = "\n".join(header)
    if NON_GRCH38_BUILD.search(text):
        raise ValueError(f"GTF 헤더가 GRCh38 이 아닌 build 를 명시한다: {gtf}\n{text}")
    status = "GRCh38" if "GRCh38" in text else "unknown"
    if status == "unknown":
        print(f"[lift][WARN] GTF 헤더에서 build 를 확인하지 못했다: {gtf} (좌표 범위 검사는 계속 적용)", file=sys.stderr)
    return {"gtf_header_build": status}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain-dir", required=True, help="<chr>.chain 이 있는 HAP 디렉터리")
    ap.add_argument("--partner-chain-dir", default="", help="반대쪽 HAP chain 디렉터리 (두 HAP 의 reference 가 같은지 검사)")
    ap.add_argument("--reference-fai", default="", help="HAP FASTA 를 만든 reference 의 .fai (chain tSize 대조)")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--out-exons", default="")
    ap.add_argument("--qc-json", default="")
    ap.add_argument("--contigs", default="", help="공백 구분. 비우면 chain 파일에서 추론")
    ap.add_argument("--reuse", action="store_true")
    a = ap.parse_args()

    cd = Path(a.chain_dir)
    contigs = a.contigs.split() if a.contigs else sorted((p.stem for p in cd.glob("chr*.chain")), key=lambda x: int(x[3:]))
    inputs = {"stage": 14, "policy": "trim", "contigs": contigs,
              "gtf": file_fingerprint(a.gtf, content=True),
              "reference_fai": optional_file(a.reference_fai, content=True),
              "chains": [file_fingerprint(cd / f"{c}.chain", content=True) for c in contigs],
              "hap_fai": [file_fingerprint(cd / f"{c}.fa.fai", content=True) for c in contigs if (cd / f"{c}.fa.fai").is_file()],
              "partner_chains": [file_fingerprint(Path(a.partner_chain_dir) / f"{c}.chain", content=True) for c in contigs] if a.partner_chain_dir else [],
              "partner_fai": [file_fingerprint(Path(a.partner_chain_dir) / f"{c}.fa.fai", content=True) for c in contigs if (Path(a.partner_chain_dir) / f"{c}.fa.fai").is_file()] if a.partner_chain_dir else [],
              "out_gtf": str(Path(a.out_gtf).resolve()),
              "out_exons": str(Path(a.out_exons).resolve()) if a.out_exons else None,
              "script": file_fingerprint(__file__, content=True)}
    old = read_cache(a.qc_json) if a.qc_json else {}
    if a.reuse and old.get("schema") == SCHEMA and cache_valid(old, inputs):
        print(f"[lift] 재사용: {a.qc_json}", file=sys.stderr)
        return 0
    ref_lengths = read_fai(Path(a.reference_fai)) if a.reference_fai else None
    maps, checks = load_chain_set(cd, contigs, "self", ref_lengths)
    if a.partner_chain_dir:
        partner, pchecks = load_chain_set(Path(a.partner_chain_dir), contigs, "partner", ref_lengths)
        checks += pchecks
        for c in contigs:
            if maps[c][1]["tSize"] != partner[c][1]["tSize"]:
                raise ValueError(f"{c}: HAP 두 개의 chain tSize 가 다르다 — 서로 다른 reference 로 만든 HAP 다")

    op = gzip.open if a.gtf.endswith(".gz") else open
    build = gtf_build_check(a.gtf, op)
    part_gtf = a.out_gtf + ".part"
    n_in_total = n_in_target = n_out = n_malformed = n_trimmed = n_deleted = n_internal = 0
    exons, out_of_range = [], []
    with op(a.gtf, "rt") as fin, open(part_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line); continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                n_malformed += 1; continue
            n_in_total += 1
            entry = maps.get(f[0])
            if entry is None:
                continue                       # chr1–22 외 contig 는 쓰지 않는다
            (ts, qs, sz), h = entry
            n_in_target += 1
            s1, e1 = int(f[3]), int(f[4])
            if e1 > h["tSize"]:
                out_of_range.append(f"{f[0]}:{s1}-{e1} > tSize {h['tSize']}"); continue
            r = map_interval(s1 - 1, e1 - 1, ts, qs, sz)
            if r is None:
                n_deleted += 1; continue
            s0, e0, aligned, boundary_trimmed = r
            if aligned < e1 - s1 + 1:
                n_internal += 1
            if boundary_trimmed:
                n_trimmed += 1
                f[8] = f[8].rstrip().rstrip(";") + '; lift_trimmed "1";'
            f[3], f[4] = str(s0 + 1), str(e0 + 1)
            fout.write("\t".join(f) + "\n"); n_out += 1
            if f[2] == "exon" and a.out_exons:
                gid = get_attr(f[8], "gene_id")
                if gid:
                    exons.append((gid, f[0], s0, e0 + 1))
    if out_of_range:
        os.remove(part_gtf)
        raise ValueError(f"GTF 좌표가 chain reference 길이를 넘는다 ({len(out_of_range)}개, 예: {out_of_range[:3]}) "
                         "— GTF 와 HAP FASTA 의 reference build 가 다르다")
    os.replace(part_gtf, a.out_gtf)
    if a.out_exons:
        with open(a.out_exons + ".part", "w") as fh:
            for gid, c, s, e in exons:
                fh.write(f"{gid}\t{c}\t{s}\t{e}\n")
        os.replace(a.out_exons + ".part", a.out_exons)
    if a.qc_json:
        qc = {"schema": SCHEMA, "on_unlifted": "trim", "input_gtf": a.gtf, "input_gtf_sha256": inputs["gtf"]["sha256"], "output_gtf": a.out_gtf,
              "output_exons": a.out_exons or None, "partner_chain_dir": a.partner_chain_dir or None,
              "reference_fai": a.reference_fai or None, "reference_checks": {**build, "chains": checks},
              "requested_contigs": contigs, "chains_loaded": len(maps),
              "features_input_total": n_in_total, "features_input_target_contigs": n_in_target,
              "features_output": n_out, "features_excluded_non_target_contigs": n_in_total - n_in_target,
              "malformed_input_lines": n_malformed, "features_trimmed": n_trimmed,
              "features_with_deleted_bases": n_internal, "features_deleted_entirely": n_deleted,
              "features_dropped_unlifted": n_deleted, "exons_output": len(exons),
              **cache_fields(inputs, [Path(a.out_gtf)] + ([Path(a.out_exons)] if a.out_exons else []))}
        tmp = Path(a.qc_json + ".part")
        tmp.write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n")
        os.replace(tmp, a.qc_json)
    print(f"[lift] feature {n_in_target:,} → {n_out:,} (경계 trim {n_trimmed:,}, 전체 결실 {n_deleted:,}) / "
          f"exon {len(exons):,}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
