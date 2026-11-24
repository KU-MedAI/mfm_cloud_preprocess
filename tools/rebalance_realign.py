#!/usr/bin/env python3
"""실행 중 BWA/HiFi minimap2의 RAM 예약을 검증함. --apply일 때만 예약을 감소함."""
from pathlib import Path
import argparse
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from integrate_wgs_rna.realign.live_budget import propose,apply


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pool-dir',type=Path,required=True)
    parser.add_argument('--cohort-root',type=Path,required=True)
    parser.add_argument('--sample',action='append',required=True)
    parser.add_argument('--apply',action='store_true')
    parser.add_argument('--report',type=Path)
    args=parser.parse_args();results=[];errors=[]
    for sample in args.sample:
        matches=[]
        for path in (args.pool_dir/'leases').glob('*.json'):
            try:row=json.loads(path.read_text())
            except (OSError,ValueError):continue
            if row.get('sample')==sample and row.get('realign_phase')=='align_restore_sort_qc':matches.append(path)
        try:
            if len(matches)!=1:raise ValueError('생존정렬lease가단일하지않음')
            plan=propose(matches[0],args.cohort_root)
            results.append(apply(plan) if args.apply else plan)
        except (OSError,ValueError,KeyError) as error:errors.append({'sample':sample,'error':str(error)})
    report={'status':'pass' if not errors else 'partial_or_blocked','apply':args.apply,'results':results,'errors':errors}
    text=json.dumps(report,ensure_ascii=False,indent=2)+'\n'
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(text)
    print(text)
    return 1 if errors else 0


if __name__=='__main__':raise SystemExit(main())
