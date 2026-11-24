#!/usr/bin/env python3
"""[구번호] 참조 좌표 GTF/exon 을 chain 으로 **개인 좌표**로 옮긴다.

현재 v2 runner에서는 이 구현을 직접 호출하지 않는다. 현재 대응은 Stage 14이며,
실행 계약은 `14_lift_annotation.sh`에서 확인한다.

    python 13_lift_annotation.py --chain-dir <individual_fasta/HAP1> --gtf <gencode.gtf.gz> \
        --out-gtf <hap1.gtf> [--out-exons <exons.hap1.tsv>]

★ 왜 필요한가
  개인 유전체에 STAR index 를 만들 때 splice junction 주석을 주려면 **좌표가 개인 기준**이어야
  한다. 참조 좌표 GTF 를 그대로 주면 indel 만큼 어긋난 junction 이 삽입돼 정렬이 나빠진다.
  ASE 의 exon 구간도 같은 이유로 lift 가 필요하다(hap1/hap2 좌표가 서로 다르다).

★ chain 해석
  `bcftools consensus --chain` 은 target=참조 → query=개인 방향으로 쓴다.
  헤더: chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
  블록 줄: size dt dq  (dt=참조 쪽 gap, dq=개인 쪽 gap)
  ⇒ 참조 좌표 t 를 개인 좌표 q 로 옮기는 구간별 offset 표를 만든다.

★ 옮길 수 없는 좌표 — `--on-unlifted` (기본 drop)
  참조 좌표가 **결실 구간 안**에 떨어지면 개인 서열에 대응 위치가 없다.

  drop (기본, AlleleSeq2 기준)
    그 feature 를 버리고 `features_dropped_unlifted` 로 센다. AlleleSeq2 는 lift 되지 않은
    site 를 `not_lifted?` 로 표시했다가 `filter_phase_warnings.py` 에서 제외한다.
    ASE 에서 스냅이 위험한 이유는 **스냅이 HAP1/HAP2 에서 비대칭으로 일어나기 때문**이다.
    두 haplotype 의 결실 위치가 다르므로 같은 exon 이 한쪽에서만 밀리고, 그러면 같은 fragment 가
    hap 별로 다른 유전자에 귀속될 수 있다 — het 와 무관한 신호가 대립 증거로 둔갑한다.

  snap (이전 동작)
    가장 가까운 블록 경계로 옮기고 `lift_snapped=1` 플래그를 남긴다. 좌표를 보존해야 하는
    비-ASE 용도(예: 커버리지 트랙 대조)에서만 쓴다.
"""
from __future__ import annotations
import argparse, bisect, gzip, json, sys
from pathlib import Path


def load_chain(path: Path):
    """(t_starts, q_starts, sizes) — 참조→개인 정렬 블록."""
    ts, qs, sz = [], [], []
    with open(path) as fh:
        hdr = fh.readline().split()
        if not hdr or hdr[0] != "chain":
            raise ValueError(f"{path}: chain 헤더가 아니다")
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
    return ts, qs, sz


def lift(pos0: int, ts, qs, sz):
    """참조 0-based → 개인 0-based. (좌표, snapped) 반환."""
    i = bisect.bisect_right(ts, pos0) - 1
    if i < 0:
        return qs[0], True
    off = pos0 - ts[i]
    if off < sz[i]:
        return qs[i] + off, False
    # 블록 밖(참조 쪽 결실 구간) → 다음 블록 시작으로 스냅
    if i + 1 < len(qs):
        return qs[i + 1], True
    return qs[i] + sz[i] - 1, True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain-dir", required=True, help="chr*.chain 이 있는 디렉토리(HAP1/HAP2)")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--out-exons", default="")
    ap.add_argument("--qc-json", default="", help="GTF lift 입력·출력·snap QC JSON 경로")
    ap.add_argument("--contigs", default="", help="공백 구분. 비우면 chain 파일에서 추론")
    ap.add_argument("--on-unlifted", choices=("drop", "snap"), default="drop",
                    help="참조 좌표가 결실 구간에 떨어질 때. drop=버린다(AlleleSeq2 기준, 기본), "
                         "snap=가장 가까운 블록 경계로 옮기고 lift_snapped=1을 남긴다")
    a = ap.parse_args()

    cd = Path(a.chain_dir)
    contigs = a.contigs.split() if a.contigs else sorted(
        (p.stem for p in cd.glob("chr*.chain")), key=lambda x: int(x[3:]))
    maps = {}
    missing_chains = []
    for c in contigs:
        p = cd / f"{c}.chain"
        if p.exists():
            maps[c] = load_chain(p)
        else:
            missing_chains.append(c)
    if missing_chains:
        raise FileNotFoundError(
            "GTF lift에 필요한 chain이 없습니다: " + ",".join(missing_chains))
    print(f"[lift] chain 로드: {len(maps)} contig", file=sys.stderr)

    op = gzip.open if a.gtf.endswith(".gz") else open
    n_in_total = n_in_target = n_out = n_snap = n_malformed = n_unlifted = 0
    exons = []
    with op(a.gtf, "rt") as fin, open(a.out_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line); continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                n_malformed += 1
                continue
            n_in_total += 1
            c = f[0]
            m = maps.get(c)
            if m is None:
                continue                       # 성염색체·scaffold 등은 버린다(autosome 전용)
            n_in_target += 1
            s0, snap1 = lift(int(f[3]) - 1, *m)
            e0, snap2 = lift(int(f[4]) - 1, *m)
            if e0 < s0:
                e0 = s0
            if snap1 or snap2:
                if a.on_unlifted == "drop":
                    # AlleleSeq2 filter_phase_warnings.py 의 not_lifted? 제외와 같은 정책.
                    n_unlifted += 1
                    continue
                n_snap += 1
                f[8] = f[8].rstrip().rstrip(";") + '; lift_snapped "1";'
            f[3], f[4] = str(s0 + 1), str(e0 + 1)
            fout.write("\t".join(f) + "\n"); n_out += 1
            if f[2] == "exon" and a.out_exons:
                gid = ""
                for kv in f[8].split(";"):
                    kv = kv.strip()
                    if kv.startswith("gene_id"):
                        gid = kv.split('"')[1]; break
                if gid:
                    exons.append((gid, c, s0, e0 + 1))

    if a.out_exons and exons:
        with open(a.out_exons, "w") as fh:
            for gid, c, s, e in exons:
                fh.write(f"{gid}\t{c}\t{s}\t{e}\n")
    snap_fraction = n_snap / n_out if n_out else 0.0
    if a.qc_json:
        qc = {
            "schema": "mfm_gtf_lift_qc_v2",
            "on_unlifted": a.on_unlifted,
            "input_gtf": a.gtf,
            "output_gtf": a.out_gtf,
            "output_exons": a.out_exons or None,
            "requested_contigs": contigs,
            "chains_loaded": len(maps),
            "features_input_total": n_in_total,
            "features_input_target_contigs": n_in_target,
            "features_output": n_out,
            "features_excluded_non_target_contigs": n_in_total - n_in_target,
            "malformed_input_lines": n_malformed,
            "features_snapped": n_snap,
            "snap_fraction_of_output": snap_fraction,
            "features_dropped_unlifted": n_unlifted,
            "unlifted_fraction_of_target": (n_unlifted / n_in_target) if n_in_target else 0.0,
            "exons_output": len(exons),
        }
        Path(a.qc_json).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n")
    print(f"[lift] feature 대상 {n_in_target:,} → {n_out:,} (전체 입력 {n_in_total:,}, "
          f"스냅 {n_snap:,} = {100*snap_fraction:.3f}%, "
          f"unlifted {a.on_unlifted} {n_unlifted:,})", file=sys.stderr)
    print(f"[lift] exon {len(exons):,} → {a.out_exons or '(생략)'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
