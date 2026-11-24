"""모델 공통 Zarr 입력·DNA 정답 마스크 계약. 원 관측 깊이는 보존함."""
from pathlib import Path
import json
import numpy as np
import zarr


def dna_target_mask(codes, actual_length):
    x=np.asarray(codes)
    if x.ndim!=1 or not 0<actual_length<=x.size or np.any(x>5) or np.any(x<0):raise ValueError('DNA code/실제 길이 오류')
    if np.any(x[:actual_length]==5) or np.any(x[actual_length:]!=5):raise ValueError('DNA padding 위치 오류')
    return (x<4)&(np.arange(x.size)<actual_length)


def masked_mean(values, mask):
    x,m=np.asarray(values),np.asarray(mask,dtype=bool)
    if x.shape!=m.shape:raise ValueError('loss/정답 mask shape 불일치')
    return float(x[m].sum(dtype=np.float64)/m.sum()) if m.any() else 0.0


def torch_dna_cross_entropy(logits, codes, observation_mask=None):
    """학습 환경용 optional adapter. 전부 무시된 batch도 미분 가능한0을 반환함."""
    import torch
    import torch.nn.functional as F
    if logits.shape[:-1]!=codes.shape or logits.shape[-1]<4:raise ValueError('logits/DNA shape 불일치')
    if bool(((codes<0)|(codes>5)).any()):raise ValueError('지원하지 않는 DNA code')
    mask=codes<4
    if observation_mask is not None:
        if observation_mask.shape!=codes.shape:raise ValueError('관측 mask shape 불일치')
        mask=mask&observation_mask.bool()
    labels=codes.long().masked_fill(~mask,-100)
    losses=F.cross_entropy(logits.reshape(-1,logits.shape[-1]),labels.reshape(-1),ignore_index=-100,reduction='none').reshape(codes.shape)
    return losses.sum()/mask.sum().clamp_min(1)


def load_window(sample_root, contig, hap, index, *, require_validated=True):
    root=Path(sample_root).resolve()
    if hap not in {'HAP1','HAP2'} or index<0:raise ValueError('HAP/window index 오류')
    if '/' in contig or contig not in {f'chr{i}' for i in range(1,23)}:raise ValueError('상염색체 이름 오류')
    if require_validated:
        qc=json.loads((root/'qc/final_package_validation.json').read_text())
        if qc.get('status')!='pass' or qc.get('failed_checks') or qc.get('sample')!=root.name:
            raise ValueError('최종 검증을 통과한 모델 입력이 아님')
    manifests=list((root/'compact').glob('*/manifest.json'))
    if len(manifests)!=1:raise ValueError('모델 layout 개수 오류')
    manifest=json.loads(manifests[0].read_text())
    if manifest.get('status')!='complete' or manifest.get('validation_level')!='full':raise ValueError('미완료 compact')
    if manifest.get('wgs_payload')=='pileup_feature':
        from .window_reader import WindowReader
        reader=WindowReader(root)
        try:return reader.load(contig,hap,index)
        finally:reader.close()
    d=manifests[0].parent/contig
    sequence=zarr.open_group(str(d/'hap_sequence.zarr'),mode='r')
    if sequence.attrs.get('sample')!=root.name:raise ValueError('DNA sample 불일치')
    g=sequence[hap]
    if index>=g['lengths'].shape[0]:raise IndexError(index)
    begin,end=int(g['offsets'][index]),int(g['offsets'][index+1])
    packed=g['bases_packed'][begin:end]
    codes=np.empty(packed.size*2,dtype=np.uint8);codes[0::2]=packed>>4;codes[1::2]=packed&15
    padded=int(g['lengths'][index]);actual=int(g['actual_length'][index]);start=int(g['raw_start'][index])
    codes=codes[:padded]
    mask=dna_target_mask(codes,actual)
    depth=zarr.open_group(str(d/'depth.zarr'),mode='r')[hap]
    lo,hi=int(depth['offsets'][index]),int(depth['offsets'][index+1])
    values=depth['depth'][lo:hi]
    if values.shape!=codes.shape or values.dtype!=np.dtype('uint32'):raise ValueError('WGS 깊이 축/형식 불일치')
    result={'dna_codes':codes,'dna_target_mask':mask,'observation_mask':np.arange(padded)<actual,
            'wgs_depth':values,'actual_length':actual,'raw_start':start,'rna_depth':None}
    labels=list((root/'rna_personal').glob('*'))
    if labels:
        if len(labels)!=1:raise ValueError('RNA label 개수 오류')
        rna=zarr.open_group(str(labels[0]/hap/'coverage_windows.zarr'),mode='r')
        if rna.attrs.get('sample')!=root.name or rna.attrs.get('hap')!=hap:raise ValueError('RNA sample/HAP 불일치')
        hits=np.flatnonzero((rna['contig'][:]==contig)&(rna['personal_start'][:]==start)&(rna['personal_end'][:]==start+actual))
        if hits.size!=1:raise ValueError('RNA/WGS 좌표 대응 오류')
        ri=int(hits[0]);rb,re=int(rna['offsets'][ri]),int(rna['offsets'][ri+1])
        observed=rna['depth'][rb:re]
        if observed.shape!=(actual,):raise ValueError('RNA 실제 길이 불일치')
        aligned=np.zeros(padded,dtype=observed.dtype);aligned[:actual]=observed
        result['rna_depth']=aligned
    return result
