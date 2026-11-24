#!/usr/bin/env python3
"""완료 패키지의 DNA/WGS/RNA 입력 대응을 읽기 전용으로 검사함."""
from pathlib import Path
import argparse
import json
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from integrate_wgs_rna.model_input import load_window


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sample_root',type=Path)
    parser.add_argument('--contig',default='chr22')
    parser.add_argument('--hap',choices=['HAP1','HAP2'],default='HAP1')
    parser.add_argument('--index',type=int,default=0,help='저장 배열에서의 0 기반 순번')
    args=parser.parse_args()
    result=load_window(args.sample_root,args.contig,args.hap,args.index)
    print(json.dumps({'status':'pass','sample':args.sample_root.name,'contig':args.contig,'hap':args.hap,
        'stored_index':args.index,'actual_length':result['actual_length'],'raw_start':result['raw_start'],
        'input_bases':len(result['dna_codes']),'dna_targets':int(result['dna_target_mask'].sum()),
        'n_count':int(np.count_nonzero(result['dna_codes']==4)),
        'padding_count':int(np.count_nonzero(result['dna_codes']==5)),
        'wgs_dtype':str(result['wgs_depth'].dtype),
        'rna_shape':list(result['rna_depth'].shape) if result['rna_depth'] is not None else None},ensure_ascii=False))


if __name__=='__main__':main()
