"""단계별 통계와 최종 모델 입력 용량을 원본 로그/QC와 함께 누적한다."""
from __future__ import annotations
import argparse, csv, fcntl, hashlib, json, os, tempfile, time
from pathlib import Path
from .cleanup import KEEP_RNA


def _tree_bytes(directory, *, rna=False):
    total=0;pending=[(Path(directory),not rna)]
    while pending:
        path,keep=pending.pop()
        if not path.is_dir():continue
        with os.scandir(path) as entries:
            for entry in entries:
                if entry.is_dir(follow_symlinks=False):
                    pending.append((Path(entry.path),keep or entry.name in KEEP_RNA))
                elif entry.is_file(follow_symlinks=False) and (keep or entry.name.endswith('.tsv')):
                    total+=entry.stat(follow_symlinks=False).st_size
    return total


def _compact_contract(root):
    manifests=sorted((root/'compact').glob('*/manifest.json'))
    if not manifests:return None
    paths=[]
    for path in manifests:
        manifest=json.loads(path.read_text())
        if manifest.get('status')!='complete':return None
        paths.append(path)
        for record in manifest.get('records',{}).values():
            directory=path.parent/str(record['contig'])
            for name in ('contig_metadata.json','filter_qc.json','contig_timing.json'):
                p=directory/name
                if p.exists():paths.append(p)
    # 완료 compact는 다음 생성 시 manifest가 running으로 바뀌는 불변 자산 계약임.
    # 임의의 chunk 변조 검사는 이 통계 cache가 아닌 기존 전수 QC가 담당함.
    digest=hashlib.sha256()
    for p in paths:
        stat=p.stat()
        digest.update(str((str(p),stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns)).encode())
        digest.update(p.read_bytes())
    return digest.hexdigest()


def _compact_bytes(root):
    metadata=root/'metadata';metadata.mkdir(parents=True,exist_ok=True)
    cache=metadata/'compact_bytes_cache.json'
    with (metadata/'compact_bytes_cache.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        try:contract=_compact_contract(root)
        except (OSError,ValueError,KeyError):contract=None
        try:saved=json.loads(cache.read_text())
        except (OSError,ValueError):saved={}
        if (contract and saved.get('schema')=='mfm_compact_bytes_cache_v1' and saved.get('contract')==contract
                and type(saved.get('bytes')) is int and saved['bytes']>=0):
            return saved['bytes'],True
        size=_tree_bytes(root/'compact')
        try:
            if contract and contract==_compact_contract(root):
                with tempfile.NamedTemporaryFile('w',dir=metadata,delete=False) as out:
                    json.dump(dict(schema='mfm_compact_bytes_cache_v1',contract=contract,bytes=size),out)
                    temporary=Path(out.name)
                try:os.replace(temporary,cache)
                finally:temporary.unlink(missing_ok=True)
        except (OSError,ValueError,KeyError):pass
        return size,False


def statistics_model_bytes(root):
    root=Path(root)
    compact,hit=_compact_bytes(root)
    return compact+_tree_bytes(root/'rna_personal',rna=True),hit

def distribution(values):
    import numpy as np
    a=np.asarray(values)
    if not a.size:return {'count':0}
    unique,counts=np.unique(a,return_counts=True)
    # quality/MAPQ는 discrete PMF, read length/depth는 histogram density를 제공한다.
    if len(unique)<=100:
        histogram=[{'value':float(v),'count':int(n),'probability':float(n/a.size)} for v,n in zip(unique,counts)]
    else:
        n,edges=np.histogram(a,bins=40)
        histogram=[{'low':float(edges[i]),'high':float(edges[i+1]),'count':int(k),
            'density':float(k/(a.size*(edges[i+1]-edges[i])))} for i,k in enumerate(n)]
    return {'count':int(a.size),'min':float(a.min()),'max':float(a.max()),'mean':float(a.mean()),
        'p01_p50_p99':[float(v) for v in np.quantile(a,[.01,.5,.99])],'histogram':histogram}

def alignment_stats(path, reference=None, limit=20000):
    import pysam
    lengths=[];mapq=[];qualities=[]
    kwargs={'reference_filename':str(reference)} if reference else {}
    with pysam.AlignmentFile(str(path),'rc' if str(path).endswith('.cram') else 'rb',**kwargs) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary:continue
            lengths.append(r.query_length or 0);mapq.append(r.mapping_quality)
            if r.query_qualities is not None and len(qualities)<2_000_000:
                qualities.extend(r.query_qualities[:2_000_000-len(qualities)])
            if len(lengths)>=limit:break
    return {'method':'파일 순서의 첫 primary read 표본; 전장 범위/무작위 추정이 아님',
      'record_limit':limit,'quality_base_limit':2_000_000,'read_length':distribution(lengths),'mapq':distribution(mapq),
      'base_quality_phred':distribution(qualities)}

def compact_stats(root):
    import numpy as np,zarr
    result=[]
    for manifest in sorted((Path(root)/'compact').glob('*/manifest.json')):
        m=json.loads(manifest.read_text())
        for rec in m.get('records',{}).values():
            contig=rec.get('contig');directory=manifest.parent/str(contig)
            row={'contig':contig,'status':rec.get('status'),'depth_summary_exact':rec.get('depth_summary'),
               'read_haplotype_counts':rec.get('read_haplotype_counts'),'shared_read_count':rec.get('shared_read_count')}
            for name in ['filter_qc','contig_metadata','contig_timing']:
                p=directory/(name+'.json')
                if p.is_file():row[name]=json.loads(p.read_text())
            depth=directory/'depth.zarr'
            if depth.is_dir():
                group=zarr.open_group(str(depth),mode='r');hist={}
                for h in ['HAP1','HAP2']:
                    if h not in group:continue
                    a=group[h]['depth'];n=a.shape[0];blocks=[]
                    for start in sorted({0,max(0,n//2-2048),max(0,n-4096)}):blocks.append(a[start:min(n,start+4096)])
                    hist[h]={'method':'window depth 배열 시작/중앙/끝 각 최대4096좌표 표본',
                        'array_positions':n,'distribution':distribution(np.concatenate(blocks))}
                row['depth_density_sampled']=hist
            result.append(row)
    return {'contigs':result,'quality_note':'compact Zarr는 Phred quality와 MAPQ를 저장하지 않는다. BAM 삭제 전 표본을 기록한다.'}

def snapshot(root, sample, stage, branch, bam=None, reference=None):
    root=Path(root);qcs=[]
    for directory,dirs,files in os.walk(root,followlinks=False):
        dirs[:]=[name for name in dirs if name!='logs' and not name.endswith('.zarr')]
        for name in files:
            if not name.endswith('.json'):continue
            p=Path(directory)/name
            if ('qc' in p.parts or name.endswith('_qc.json')) and p.stat().st_size<=2_000_000:
                qcs.append(str(p.relative_to(root)))
    d={'sample':sample,'stage':stage,'branch':branch,'time':time.strftime('%Y-%m-%dT%H:%M:%S%z'),
        'qc_files':sorted(qcs),'measurement_scope':'아직 생성되지 않은 통계는 누락 상태로 둔다.'}
    if bam:
        try:d['alignment_sample']=alignment_stats(bam,reference)
        except Exception as e:d['alignment_sample']={'status':'missing','reason':str(e)}
    if stage in {'10','22','23','24','pilot'}:
        d['model_input_bytes'],hit=statistics_model_bytes(root)
        d['model_bytes_measurement']={'completed_compact_cache_hit':hit,'rna_files':'current_stat',
            'scope':'완료 compact 불변 manifest 계약의 통계 cache임. QC/cleanup 전수 집계와 독립임'}
        try:d['compact']=compact_stats(root)
        except Exception as e:d['compact']={'status':'missing','reason':str(e)}
    target=root/'metadata/stage_statistics.jsonl';target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('a') as f:
        fcntl.flock(f,fcntl.LOCK_EX);f.write(json.dumps(d,ensure_ascii=False)+'\n')
    return d

if __name__=='__main__':
    p=argparse.ArgumentParser(description='단계별 통계 누적')
    for name in ['root','sample','stage','branch']:p.add_argument('--'+name,required=True)
    p.add_argument('--bam');p.add_argument('--reference');a=p.parse_args()
    snapshot(a.root,a.sample,a.stage,a.branch,a.bam,a.reference)
