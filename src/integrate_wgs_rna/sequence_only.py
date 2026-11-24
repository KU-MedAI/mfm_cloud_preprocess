"""기존 window 원장을 그대로 사용하며 WGS 서열만 Zarr로 저장함."""
from __future__ import annotations

import argparse,csv,hashlib,json,multiprocessing,os,shutil,time,uuid
from pathlib import Path

import numpy as np
import pysam,zarr
from numcodecs import Blosc,blosc
from pre_training.preprocess_raw_readset.windowing import sequence_annotations

from .cpu_executor import cpu_results
from .legacy_shared_hap_converter import SEQUENCE_ALLOWED,pack_4bit_acgtnp_chunk
from .pipeline.checkpoints import atomic_json,digest,file_identity
from .validation_parallel import sequence_checks
from .window_contract import read_windows,validate_n_values
from .window_policy import persist_window_annotations,validate_window_annotations

MODE='sequence_only'
LAYOUT='sequence_4bit_acgtnp__windowed'
# 기존 출력의 생산 지문을 유지함. 서열 생산 규격 변경 시에만 버전을 갱신함.
PRODUCER_REVISION = '6a7950837ecdc106aed6fca29543c2376658ac75bca2d29eb723246228887a12'
PACKED_CHUNK=32<<20
PARENT_BYTES=1<<30
WORKER_BYTES=2<<30
TABLE=np.frombuffer(b'ACGTNP',dtype=np.uint8)


def _hash_window(state,start,actual,sequence):
    state.update(f'{start}:{actual}:{len(sequence)}\n'.encode('ascii'))
    state.update(sequence.encode('ascii'))


def _decode(array,begin,end,length):
    packed=np.asarray(array[begin:end]);codes=np.empty(len(packed)*2,dtype=np.uint8)
    codes[::2]=packed>>4;codes[1::2]=packed&15
    if np.any(codes[:length]>5):raise ValueError('서열4bit 코드 범위 오류')
    return TABLE[codes[:length]].tobytes().decode('ascii')


def validate_contig(directory,*,source=False,payload=MODE):
    """저장 배열·좌표·내용 SHA를 전수검사하며 source=True는 FASTA와 독립 대조함."""
    directory=Path(directory);meta=json.loads((directory/'contig_metadata.json').read_text())
    if meta.get('wgs_payload')!=payload:raise ValueError('서열 payload 계약 오류')
    if payload==MODE and (directory/'pileup.zarr').exists():raise ValueError('서열 전용 입력에 pileup이 섞임')
    if any((directory/n).exists() for n in ['depth.zarr','readset.zarr','readset.HAP1.zarr','readset.HAP2.zarr']):
        raise ValueError('서열 전용 입력에 WGS read/depth 산출물이 섞임')
    group=zarr.open_group(str(directory/'hap_sequence.zarr'),mode='r')
    if group.attrs.get('sample')!=meta['sample_id'] or group.attrs.get('contig')!=meta['contig']:
        raise ValueError('서열 sample/contig 불일치')
    rows=read_windows(directory/'windows.tsv',compact=True);result={}
    for hap in ['HAP1','HAP2']:
        chosen=[r for r in rows if r['haplotype']==hap];g=group[hap]
        if set(g.array_keys())!={'bases_packed','offsets','lengths','actual_length','raw_start'}:raise ValueError('서열 배열 계약 오류')
        offsets=g['offsets'][:];lengths=g['lengths'][:];starts=g['raw_start'][:];actual=g['actual_length'][:]
        n=len(chosen)
        if offsets.shape!=(n+1,) or any(a.shape!=(n,) for a in [lengths,starts,actual]):raise ValueError('서열 좌표 shape 불일치')
        if offsets[0]!=0 or offsets[-1]!=g['bases_packed'].shape[0] or np.any(np.diff(offsets)!=(lengths.astype(np.int64)+1)//2):
            raise ValueError('서열 packed offset 불일치')
        observed=hashlib.sha256();fasta=None
        if source:fasta=pysam.FastaFile(meta['sources'][hap]['fasta'])
        try:
            for i,row in enumerate(chosen):
                if int(starts[i])!=int(row['raw_start']) or int(actual[i])!=int(row['actual_length']) or not 0<actual[i]<=lengths[i]:
                    raise ValueError('서열 window 좌표 불일치')
                seq=_decode(g['bases_packed'],int(offsets[i]),int(offsets[i+1]),int(lengths[i]))
                if seq[int(actual[i]):]!='P'*(int(lengths[i])-int(actual[i])):raise ValueError('window padding 불일치')
                bases=seq[:int(actual[i])]
                validate_n_values(row)
                annotations=sequence_annotations(bases)
                if any(annotations[k]!=int(row[k]) for k in ['n_count','n_nonacgt','n_max_run']):
                    raise ValueError('서열과 window N 주석 불일치')
                if fasta and bases!=fasta.fetch(meta['contig'],int(starts[i]),int(starts[i]+actual[i])).upper():
                    raise ValueError('원본 FASTA와 서열 Zarr 불일치')
                _hash_window(observed,int(starts[i]),int(actual[i]),seq)
        finally:
            if fasta:fasta.close()
        if observed.hexdigest()!=meta['haplotypes'][hap]['window_content_sha256']:raise ValueError('서열 내용 SHA 불일치')
        result[hap]=dict(status='pass',level='full',windows=n,window_content_sha256=observed.hexdigest())
    return result


def _contig_task(root,layout,contig,rows,window_size,payload=MODE):
    """한 contig의2HAP을 bounded 배치로 작성하며 최종 공개는 원자적으로 수행함."""
    blosc.set_nthreads(1);root=Path(root);layout=Path(layout);destination=layout/contig
    sources={};identities=[]
    for hap in ['HAP1','HAP2']:
        fasta=root/'individual_fasta'/hap/(contig+'.fa');chain=fasta.with_suffix('.chain')
        sources[hap]=dict(fasta=str(fasta),chain=str(chain) if chain.is_file() else '')
        identities.extend(file_identity(p) for p in [fasta,Path(str(fasta)+'.fai')]+([chain] if chain.is_file() else []))
    fingerprint=digest(dict(inputs=identities,rows=rows,window_size=window_size,producer=PRODUCER_REVISION))
    if destination.exists():
        meta=json.loads((destination/'contig_metadata.json').read_text())
        if meta['input_provenance']['fingerprint_sha256']!=fingerprint:raise ValueError('기존 서열과 입력 지문이 다름; 기존 결과를 덮어쓰지 않음')
        validate_contig(destination,source=True,payload=payload)
        return dict(status='pass',sample=root.name,contig=contig,zarr=str(destination),reused=True)
    temporary=layout/(contig+'.partial.'+uuid.uuid4().hex);temporary.mkdir()
    try:
        codec=Blosc(cname='zstd',clevel=5,shuffle=Blosc.BITSHUFFLE)
        g=zarr.open_group(str(temporary/'hap_sequence.zarr'),mode='w')
        g.attrs.update(schema='mfm_hap_sequence_4bit_acgtnp_direct_v0',sample=root.name,contig=contig,
                       window_size=window_size,encoding='4bit_acgtnp',symbols='A,C,G,T,N,P',wgs_payload=payload)
        compact_rows=[];hap_meta={}
        for hap in ['HAP1','HAP2']:
            selected=[r for r in rows if r['haplotype']==hap];n=len(selected);h=g.create_group(hap)
            h.attrs.update(encoding='4bit_acgtnp',kind='hap_window_sequence',bases_per_byte=2,coordinate_system='individual_raw_0based_halfopen')
            lengths=np.full(n,window_size,dtype=np.int32);offsets=np.arange(n+1,dtype=np.int64)*((window_size+1)//2)
            values=dict(lengths=lengths,offsets=offsets,raw_start=np.array([int(r['raw_start']) for r in selected],dtype=np.int64),
                        actual_length=np.array([int(r['actual_length']) for r in selected],dtype=np.int64))
            for name,value in values.items():h.create_dataset(name,data=value,chunks=(1024,),compressor=codec)
            array=h.create_dataset('bases_packed',shape=(int(offsets[-1]),),chunks=(PACKED_CHUNK,),dtype='u1',compressor=codec)
            content=hashlib.sha256();batch=max(1,(64<<20)//window_size)
            with pysam.FastaFile(sources[hap]['fasta']) as fasta:
                for first in range(0,n,batch):
                    selected_batch=selected[first:first+batch];sequences=[]
                    for i,row in enumerate(selected_batch,first):
                        start=int(row['raw_start']);actual=int(row['actual_length'])
                        if actual>window_size:raise ValueError('window 크기보다 실제 길이가 큼')
                        seq=fasta.fetch(contig,start,start+actual).upper()
                        if len(seq)!=actual:raise ValueError('FASTA window 범위 오류')
                        seq+='P'*(window_size-actual);sequences.append(seq);_hash_window(content,start,actual,seq)
                        compact_rows.append(dict(row,window_index=i,sequence_offset=int(offsets[i]),sequence_packed_bytes=int(offsets[i+1]-offsets[i])))
                    packed,_,_,_,invalid=pack_4bit_acgtnp_chunk(sequences,SEQUENCE_ALLOWED,first,'fail')
                    if invalid:raise ValueError('FASTA에 허용되지 않은 염기 존재')
                    array[int(offsets[first]):int(offsets[first])+len(packed)]=packed
            hap_meta[hap]=dict(window_count=n,window_content_sha256=content.hexdigest())
        if not compact_rows:raise ValueError('contig에 저장할 window가 없음')
        with (temporary/'windows.tsv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(compact_rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(compact_rows)
        meta=dict(sample_id=root.name,contig=contig,wgs_payload=payload,validation_level='full',sources=sources,haplotypes=hap_meta,
                  input_provenance=dict(contig=contig,validation_level='full',fingerprint_sha256=fingerprint))
        atomic_json(temporary/'contig_metadata.json',meta)
        meta['validation']=dict(level='full',by_haplotype=validate_contig(temporary,payload=payload))
        atomic_json(temporary/'contig_metadata.json',meta)
        if identities!=[file_identity(item[0]) for item in identities]:raise ValueError('서열 작성 중 FASTA/chain 입력 변경')
        atomic_json(temporary/'filter_qc.json',dict(wgs_payload=payload,wgs_read_filters='not_applicable' if payload==MODE else 'pileup_policy',window_index_selection='unchanged',kept_windows=len(compact_rows)))
        os.rename(temporary,destination)
    except BaseException:
        shutil.rmtree(temporary,ignore_errors=True);raise
    return dict(status='pass',sample=root.name,contig=contig,zarr=str(destination))


def export_sequence(root,*,contigs,jobs=None):
    root=Path(root).resolve();index=root/'windows';rows=read_windows(index/'windows.tsv')
    config=json.loads((index/'windows.json').read_text());window_size=int(config['window_bp'])
    selected={c:[r for r in rows if r['contig']==c and r['kept'].lower() in {'1','true'}] for c in contigs}
    if set(contigs)!={r['contig'] for r in rows}:raise ValueError('전체 contig/window 계약 불일치')
    layout=root/'compact'/LAYOUT;layout.mkdir(parents=True,exist_ok=True)
    records={};manifest=dict(schema='mfm_sequence_only_manifest_v1',sample=root.name,wgs_payload=MODE,validation_level='full',status='running',records=records)
    manifest['memory_plan']={'parent_bytes':PARENT_BYTES,'worker_bytes':WORKER_BYTES,'batch_base_bytes':64<<20,'worker_address_space_bytes':WORKER_BYTES}
    atomic_json(layout/'manifest.json',manifest)
    cpu=max(1,int(os.environ.get('MFM_SAMPLE_CPU_BUDGET','1')));upper=max(1,int(os.environ.get('MFM_CPU_BUDGET',cpu)))
    if jobs is not None:upper=min(upper,max(1,int(jobs)))
    count=min(len(contigs),max(1,upper-1));old=os.environ.get('MFM_MEMORY_STAGE_BASE_BYTES');os.environ['MFM_MEMORY_STAGE_BASE_BYTES']=str(PARENT_BYTES)
    try:
        for _,result,error in cpu_results(_contig_task,[(str(root),str(layout),c,selected[c],window_size) for c in contigs],min(cpu,count),
                elastic=bool(os.environ.get('MFM_CPU_POOL_LEASE_FILE')),max_workers=count,mp_context=multiprocessing.get_context('spawn'),
                memory_per_worker_bytes=WORKER_BYTES,parent_cpus=1):
            if error:raise error
            records[result['contig']]=result;atomic_json(layout/'manifest.json',manifest)
    finally:
        if old is None:os.environ.pop('MFM_MEMORY_STAGE_BASE_BYTES',None)
        else:os.environ['MFM_MEMORY_STAGE_BASE_BYTES']=old
    annotations=[]
    for contig in contigs:
        row=dict(contig=contig)
        for hap in ['hap1','hap2']:
            p=root/'individual_fasta'/hap.upper()/(contig+'.chain');row[hap+'_chain']=str(p) if p.exists() else ''
        annotations.append(row)
    persist_window_annotations(root,layout,index,annotations,None)
    validate_window_annotations(root)
    checks = sequence_checks((layout / contig for contig in contigs),
                             source=True, max_workers=count)
    trace = {}
    for contig, (_, error) in zip(contigs, checks, strict=True):
        if error is not None:
            raise ValueError(f'{contig} 원본 FASTA 대조 실패: {error}')
        trace[contig] = 'pass'
    manifest['status']='complete';atomic_json(layout/'manifest.json',manifest)
    atomic_json(root/'qc/stage12_trace/sequence_only/trace_manifest.json',dict(schema='mfm_stage70_trace_manifest_v1',status='pass',
                sample=root.name,compact_root=str(layout),wgs_payload=MODE,contigs=trace,scope='원본FASTA와 모든kept window 서열·좌표·padding 대조',completed_at=time.time()))
    return manifest


def preflight(root):
    """read QC를 위장하지 않고 실제 개인FASTA/VCF 자산 검증을 별도로 기록함."""
    from tools.personal_assets import verify
    root=Path(root);qc=json.loads((root/'individual_fasta/personal_fasta_qc.json').read_text())
    if qc.get('schema')!='mfm_personal_fasta_sample_qc_v1' or qc.get('sample')!=root.name or qc.get('status')!='pass':
        raise ValueError('개인FASTA 생성 QC 계약 불일치')
    result=verify(root/'individual_fasta')
    if result['status']!='pass':raise ValueError('개인FASTA/chain 자산 무결성 실패')
    atomic_json(root/'qc/sequence_preflight.json',dict(schema='mfm_sequence_preflight_v1',sample=root.name,status='pass',
                scope='개인FASTA 생성QC 및 FASTA/FAI/chain SHA 검증; WGS read QC와 구분함',assets=result))


def stamp_metadata(root):
    root=Path(root);row={}
    if os.environ.get('PIPELINE_MANIFEST'):
        with Path(os.environ['PIPELINE_MANIFEST']).open() as handle:
            rows=[r for r in csv.DictReader(handle,delimiter='\t') if r['sample_id']==root.name]
        if len(rows)!=1:raise ValueError('서열 sample manifest 불일치')
        row=rows[0]
        shutil.copyfile(os.environ['PIPELINE_MANIFEST'],root/'metadata/source_manifest.tsv')
    value=dict(row,sample_id=root.name,sample=root.name,wgs_payload=MODE,has_rna=bool(row.get('rna_path')),
               wgs_readset_available=False,wgs_depth_available=False,compact_root=str(root/'compact'/LAYOUT),
               rna_payload=os.environ.get('MFM_RNA_OUTPUT_MODE','full'))
    if os.environ.get('MFM_SEQUENCE_DEFER_RNA')=='1' and value['has_rna']:
        value.update(completion_scope='dna_sequence_only_rna_deferred',rna_processing_status='deferred',
                     wgs_reads_processed=False,rna_coverage_available=False,ase_available=False,cnv_available=False)
    atomic_json(root/'metadata/sample.json',value)
    atomic_json(root/'metadata/preprocess.json',dict(schema='mfm_sequence_only_preprocess_v1',sample=root.name,wgs_payload=MODE,
                window_index=str(root/'metadata/wgs_windows.tsv'),compact_root=str(root/'compact'/LAYOUT),encoding='4bit_acgtnp',compressor='zstd',compression_level=5))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);parser.add_argument('--contigs',required=True)
    args=parser.parse_args();preflight(args.root)
    result=export_sequence(args.root,contigs=args.contigs.replace(',',' ').split());stamp_metadata(args.root)
    print(json.dumps(result,ensure_ascii=False))
