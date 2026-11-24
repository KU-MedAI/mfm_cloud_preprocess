#!/usr/bin/env python3
"""기존 완료본의 WGS window 관측 깊이를 집계하여 지정한 TSV에 저장함."""
from pathlib import Path
import argparse
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from integrate_wgs_rna.window_quality import quality_from_existing,write_quality


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('compact_root',type=Path,help='chr1~22가 있는 layout 폴더')
    parser.add_argument('--contigs',nargs='+',help='선택할 염색체; 생략하면 전체')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.resolve().is_relative_to(args.compact_root.resolve()):
        parser.error('기존 compact를 보존하도록 출력은 compact 바깥 경로로 지정해야 함')
    if args.output.exists():parser.error('기존 집계 파일을 덮어쓰지 않음; 새 출력 경로를 지정해야 함')
    write_quality(args.output,quality_from_existing(args.compact_root,set(args.contigs) if args.contigs else None))
    print(f'관측 깊이 집계 저장: {args.output}')


if __name__=='__main__':main()
