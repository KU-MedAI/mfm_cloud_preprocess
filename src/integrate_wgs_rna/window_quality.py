"""padding을 제외한 관측 깊이 원값 통계. 저장/학습 선택에는 적용하지 않음."""
from pathlib import Path
import csv
import numpy as np


def depth_metrics(values, actual_length, *, backend=None):
    import os
    if (backend or os.environ.get('MFM_QC_BACKEND','cpu'))=='cupy':
        from .gpu_arrays import depth_metrics as gpu_depth_metrics
        return gpu_depth_metrics(values,actual_length)
    values=np.asarray(values)
    if actual_length <= 0 or actual_length > values.size:raise ValueError('관측 구간 길이 오류')
    actual=values[:actual_length]
    if actual.dtype.kind not in 'ui' or np.any(actual<0):raise ValueError('음수 관측 깊이')
    total=int(actual.sum(dtype=np.uint64))
    nonzero=int(np.count_nonzero(actual))
    return {'depth_sum':total,'depth_mean':total/actual_length,
            'depth_zero_fraction':(actual_length-nonzero)/actual_length,'depth_max':int(actual.max())}


def write_quality(path, rows):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    columns=['contig','haplotype','window_index','raw_start','actual_length',
             'depth_sum','depth_mean','depth_zero_fraction','depth_max']
    with p.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns,delimiter='\t')
        writer.writeheader();writer.writerows(rows)


def quality_from_existing(compact_root, contigs=None):
    """구완료본 분석용. 배열을 window단위로 읽고 Zarr 자체는 수정하지 않음."""
    import zarr
    root=Path(compact_root)
    for d in sorted(root.glob('chr*')):
        if not d.is_dir() or (contigs and d.name not in contigs):continue
        with (d/'windows.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        group=zarr.open_group(str(d/'depth.zarr'),mode='r')
        for hap in ['HAP1','HAP2']:
            selected=[r for r in rows if r['haplotype']==hap]
            offsets=group[hap]['offsets'][:]
            for i,r in enumerate(selected):
                length=int(r['actual_length'])
                values=group[hap]['depth'][int(offsets[i]):int(offsets[i])+length]
                yield {**{k:r[k] for k in ['contig','haplotype','window_index','raw_start','actual_length']},**depth_metrics(values,length)}
