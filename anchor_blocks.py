#!/usr/bin/env python3
"""리드 기반 phase 블록을 **패널 phase 방향에 앵커링**한다.

    python anchor_blocks.py --panel <panel.vcf.gz> --read <whatshap_phased.vcf.gz> \
        --sample S --out <anchored.vcf.gz> [--stats anchored.stats.tsv]

★ 왜 필요한가
  `whatshap phase` 는 **리드가 잇는 구간(블록)** 안에서만 phase 를 정한다. 블록 사이의 상대적
  방향은 리드로 알 수 없어 임의다. 그래서 리드 phase 만으로 H1/H2 FASTA 를 만들면
  **블록 경계마다 haplotype 이 뒤바뀐 서열**이 나온다.

  해결: 블록 내부 순서는 리드를 믿고, **블록의 전역 방향은 패널로 정한다.**
  블록의 het 사이트에서 리드 phase 와 패널 phase 의 일치/불일치를 세어 다수결로 방향을 잡고,
  불일치가 많으면 블록을 통째로 뒤집는다. ⇒ H1 은 염색체 끝까지 같은 H1 이다.

  이 메커니즘은 phASER 가 `blockGWPhase`/`gwStat` 로 하는 것과 동일하다(실측 검증됨:
  gw_phased=1 이 85,485 / =0 이 347 = 99.6%).

★ 우선순위: **리드 1순위 / 패널 2순위** (jhleee 지시 2026-08-26)
  - 블록 **내부** phase 는 언제나 리드를 따른다. 패널과 다르면 리드가 이긴다 — 그게 교정의 정의다.
  - 패널은 두 가지에만 쓴다: ① 블록의 **전역 방향**(어느 쪽이 H1인가) ② 리드가 phase 하지
    못한 사이트의 **대체값**.
  ⚠ 초기 설계는 "일치율 < 0.90 이면 리드를 버린다" 였고 이는 목적과 반대로 작동했다.
    실측(HG00149 chr20): 24 Mb 블록의 일치율이 **0.53** 이다. 그 안에 switch 가 276번 있으면
    일치/불일치가 반반이 되기 때문이다. 고정 임계로는 "리드가 많이 교정한 블록" 과
    "방향을 못 정하는 블록" 을 구분할 수 없다.

★ 방향 결정은 이항 유의성으로 한다
  방향은 다수결이되, 그 다수결이 **동전 던지기와 구분되는가**를 본다:
      z = |n_agree - n_disagree| / sqrt(n_agree + n_disagree)
  z >= --orient-min-z (기본 3) 이면 방향을 확정한다. 위 0.53 블록도 |13269-11764|=1505,
  sqrt(25033)=158 → z=9.5 로 압도적이다. z 가 낮은 블록(대개 het 이 몇 개뿐)만 방향을
  못 정하고, 그 경우 뒤집지 않는다(패널 방향 유지) + 플래그를 남긴다.

출력: 패널과 같은 사이트 집합, GT 는 블록 내부에서 리드가 교정한 phase, 방향은 패널 기준.
      `PS` 는 리드 블록 id 를 보존하고, 앵커링 결정은 stats TSV 에 남긴다.
"""
from __future__ import annotations
import argparse, sys
from collections import defaultdict

import pysam


def gt_of(rec, idx=0):
    s = rec.samples[idx]
    gt = s.get("GT")
    if gt is None or None in gt or len(gt) != 2:
        return None
    return tuple(gt)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True, help="패널 phased VCF (염색체 전체 일관)")
    ap.add_argument("--read", required=True, help="whatshap phase 산출 VCF (블록 단위)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stats", default="")
    ap.add_argument("--min-conf", type=float, default=0.0,
                    help="(폐기) 리드 1순위 정책에서는 쓰지 않는다. 0 = 게이트 없음")
    ap.add_argument("--min-het", type=int, default=1)
    ap.add_argument("--orient-min-z", type=float, default=3.0,
                    help="블록 방향 확정에 필요한 이항 z. 미달이면 뒤집지 않고 플래그만 남긴다")
    ap.add_argument("--mode", default="anchored", choices=["anchored", "readonly"],
                    help="anchored=패널 골격+리드 정제 | readonly=리드 phase 만(앵커링 없음)")
    ap.add_argument("--block-table", default="", help="readonly 모드에서 블록 경계 표")
    a = ap.parse_args()

    # ── 1) 리드 phase 를 블록별로 모은다 ──
    #  whatshap 은 phased het 에 PS(phase set) 를 붙인다. PS 가 없으면 그 사이트는 unphased.
    read_gt: dict[tuple[str, int], tuple] = {}
    read_ps: dict[tuple[str, int], int] = {}
    with pysam.VariantFile(a.read) as vf:
        si = list(vf.header.samples).index(a.sample) if a.sample in vf.header.samples else 0
        for rec in vf:
            gt = gt_of(rec, si)
            if gt is None or gt[0] == gt[1]:
                continue
            s = rec.samples[si]
            if not s.phased:
                continue
            ps = s.get("PS")
            if ps is None:
                continue
            key = (rec.contig, rec.pos)
            read_gt[key] = gt
            read_ps[key] = int(ps)

    # ── 2) 블록별로 패널과의 일치/불일치 집계 ──
    agree = defaultdict(int); disagree = defaultdict(int)
    panel_gt: dict[tuple[str, int], tuple] = {}
    with pysam.VariantFile(a.panel) as vf:
        si = list(vf.header.samples).index(a.sample) if a.sample in vf.header.samples else 0
        for rec in vf:
            gt = gt_of(rec, si)
            if gt is None or gt[0] == gt[1]:
                continue
            key = (rec.contig, rec.pos)
            panel_gt[key] = gt
            rg = read_gt.get(key)
            if rg is None:
                continue
            b = read_ps[key]
            if rg == gt:
                agree[b] += 1
            else:
                disagree[b] += 1

    # ── 3) 블록 방향 결정 ──
    #  flip=True  → 리드 phase 를 뒤집어 패널 방향에 맞춘다
    #  use_read=False → 신뢰도 미달. 교정하지 않고 패널 phase 를 그대로 쓴다
    #  리드 1순위: use_read 는 **항상 True** 다(블록 내부는 리드를 따른다).
    #  판정하는 것은 "뒤집을지" 뿐이고, 그 판정이 유의하지 않으면 뒤집지 않는다.
    decision: dict[int, tuple[bool, bool, int, int, float]] = {}
    n_unoriented = 0
    for b in set(agree) | set(disagree):
        na, nd = agree[b], disagree[b]
        n = na + nd
        z = abs(na - nd) / (n ** 0.5) if n else 0.0
        oriented = (n >= a.min_het) and (z >= a.orient_min_z)
        flip = (nd > na) and oriented
        if not oriented:
            n_unoriented += 1
        decision[b] = (flip, True, na, nd, z)      # use_read=True 고정

    # ── 4a) readonly 모드: 리드가 phase 한 사이트만, 앵커링 없이 그대로 ──
    #  블록 간 방향이 임의이므로 하류가 경계를 알아야 한다 → block_table 을 함께 낸다.
    #  리드가 phase 하지 못한 het 은 **버린다**(임의 배정보다 누락이 낫다).
    if a.mode == "readonly":
        n_kept = n_drop_het = n_hom = 0
        blocks_seen: dict[int, list] = {}
        with pysam.VariantFile(a.panel) as vin:
            hdr = vin.header.copy()
            if "PS" not in hdr.formats:
                hdr.add_line('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set (read-based block)">')
            hdr.add_line('##phase_mode=readonly (no panel anchoring; H1/H2 identity is per-block)')
            with pysam.VariantFile(a.out, "w", header=hdr) as vout:
                si = list(vin.header.samples).index(a.sample) if a.sample in vin.header.samples else 0
                for rec in vin:
                    gt = gt_of(rec, si)
                    if gt is None:
                        continue
                    key = (rec.contig, rec.pos)
                    if gt[0] == gt[1]:
                        h = rec.copy(); h.translate(hdr)
                        vout.write(h); n_hom += 1; continue        # hom-alt 는 phase 무관
                    rg = read_gt.get(key)
                    if rg is None:
                        n_drop_het += 1; continue                   # 리드가 못 잇는 het → 제외
                    out = rec.copy()
                    #  ⚠ rec 은 **입력 헤더**에 묶여 있다. 그대로 PS 를 쓰면
                    #    KeyError: unknown format: PS 로 죽는다 → 출력 헤더로 재바인딩한다.
                    out.translate(hdr)
                    out.samples[0]["GT"] = rg
                    out.samples[0].phased = True
                    b = read_ps[key]; out.samples[0]["PS"] = b
                    vout.write(out); n_kept += 1
                    blk = blocks_seen.setdefault(b, [rec.contig, rec.pos, rec.pos, 0])
                    blk[2] = rec.pos; blk[3] += 1
        tot_het = len(panel_gt)
        print(f"[readonly] het {tot_het:,} 중 리드 phase {n_kept:,} 유지 / {n_drop_het:,} 제외 "
              f"({100*n_drop_het/max(tot_het,1):.1f}% 손실) / hom-alt {n_hom:,} 보존")
        print(f"[readonly] 블록 {len(blocks_seen):,}개 — 블록 경계마다 H1/H2 정체가 뒤바뀔 수 있다")
        if a.block_table:
            with open(a.block_table, "w") as fh:
                fh.write("block_ps\tcontig\tstart\tend\tn_het\tspan_bp\n")
                for b, (c, s0, e0, n) in sorted(blocks_seen.items(), key=lambda x: (x[1][0], x[1][1])):
                    fh.write(f"{b}\t{c}\t{s0}\t{e0}\t{n}\t{e0-s0+1}\n")
            spans = sorted((e0-s0+1 for _, (_, s0, e0, _) in blocks_seen.items()), reverse=True)
            if spans:
                tot = sum(spans); acc = 0
                n50 = next(x for x in spans if (acc := acc + x) >= tot/2)
                print(f"[readonly] 블록 span 합 {tot/1e6:.1f} Mb / N50 {n50/1e3:.1f} kb / 최대 {spans[0]/1e6:.2f} Mb")
            print(f"[readonly] 블록 표 → {a.block_table}")
        return 0

    # ── 4b) anchored 모드: 패널 사이트 집합 + 교정된 GT ──
    n_out = n_corrected = n_flipped_block = n_panel_kept = 0
    with pysam.VariantFile(a.panel) as vin:
        hdr = vin.header.copy()
        if "PS" not in hdr.formats:
            hdr.add_line('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set (read-based block)">')
        hdr.add_line(f'##anchor_blocks=min_conf={a.min_conf},min_het={a.min_het}')
        with pysam.VariantFile(a.out, "w", header=hdr) as vout:
            si = list(vin.header.samples).index(a.sample) if a.sample in vin.header.samples else 0
            for rec in vin:
                out = rec.copy()
                gt = gt_of(rec, si)
                if gt is not None and gt[0] != gt[1]:
                    key = (rec.contig, rec.pos)
                    rg = read_gt.get(key)
                    if rg is not None:
                        b = read_ps[key]
                        flip, use_read, *_ = decision[b]
                        if use_read:
                            out.translate(hdr)      # 출력 헤더 재바인딩(PS 쓰기 위해)
                            new = (rg[1], rg[0]) if flip else rg
                            if new != gt:
                                n_corrected += 1
                            out.samples[0]["GT"] = new
                            out.samples[0].phased = True
                            out.samples[0]["PS"] = b
                        else:
                            n_panel_kept += 1
                vout.write(out)
                n_out += 1
    n_flipped_block = sum(1 for v in decision.values() if v[0])

    tot_het = len(panel_gt)
    covered = sum(1 for k in panel_gt if k in read_gt)
    print(f"[anchor] 정책: **리드 1순위 / 패널 2순위** (방향만 패널, 내부는 리드)")
    print(f"[anchor] 패널 het {tot_het:,} / 리드가 phase 한 het {covered:,} ({100*covered/max(tot_het,1):.1f}%)")
    print(f"[anchor] 리드가 못 잇는 het {tot_het-covered:,} → 패널 phase 로 대체(2순위)")
    print(f"[anchor] 블록 {len(decision):,}개 / 뒤집은 블록 {n_flipped_block:,} / 방향 미확정(z<{a.orient_min_z}) {n_unoriented:,}")
    print(f"[anchor] ★ 최종 phase 가 패널과 다른 het: {n_corrected:,}  ({100*n_corrected/max(covered,1):.2f}% of 리드 phase 된 het)")

    if a.stats:
        with open(a.stats, "w") as fh:
            fh.write("block_ps\tn_agree\tn_disagree\tz\tflipped\toriented\n")
            for b, (flip, use_read, na, nd, z) in sorted(decision.items()):
                fh.write(f"{b}\t{na}\t{nd}\t{z:.2f}\t{int(flip)}\t{int(z >= a.orient_min_z)}\n")
        print(f"[anchor] 블록별 결정 → {a.stats}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
