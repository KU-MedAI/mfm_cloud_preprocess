#!/usr/bin/env python3
"""최종 PASS 패키지의 N/P·DNA/depth·RNA 대응을 읽기 전용으로 검증함."""
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from integrate_wgs_rna.input_audit import audit_sample


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sample_root',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int)
    parser.add_argument('--max-per-axis',type=int,default=0,help='0은 전수, 양수는 축별 표본 검사임')
    args=parser.parse_args()
    report=audit_sample(args.sample_root,args.output,workers=args.workers,max_per_axis=args.max_per_axis)
    print(json.dumps({k:report[k] for k in ['status','sample','scope','totals','wall_seconds']},ensure_ascii=False))
    return report['status']!='pass'


if __name__=='__main__':sys.exit(main())
