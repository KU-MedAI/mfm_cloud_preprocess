"""완성 BAM의 단계 QC 통과 후 중복 염색체 조각만 정리함."""
from pathlib import Path
import hashlib,json,os,re,subprocess,tempfile,time
import pysam


def fingerprint(path):
    value=path.stat()
    return [value.st_dev,value.st_ino,value.st_size,value.st_mtime_ns,value.st_ctime_ns]


def opened_files(paths):
    wanted={tuple(fingerprint(p)[:2]):str(p) for p in paths};found=[]
    for process in Path('/proc').iterdir():
        if not process.name.isdigit():continue
        try:
            for fd in (process/'fd').iterdir():
                try:
                    value=fd.stat();key=(value.st_dev,value.st_ino)
                    if key in wanted:found.append({'pid':int(process.name),'path':wanted[key]})
                except OSError:pass
        except OSError:pass
    return found


def atomic(path,value):
    fd,name=tempfile.mkstemp(prefix='.'+path.name,dir=path.parent)
    try:
        with os.fdopen(fd,'w') as handle:json.dump(value,handle,ensure_ascii=False,indent=2);handle.write('\n')
        os.replace(name,path)
    finally:Path(name).unlink(missing_ok=True)


def cleanup_parallel_shards(qc_json,samtools='samtools',*,dry_run=False):
    qc_json=Path(qc_json);qc=json.loads(qc_json.read_text());folder=qc_json.parent
    if qc.get('schema')!='mfm_haplotag_sample_qc_v1' or qc.get('status')!='pass':
        raise ValueError('하플로태깅 단계 QC PASS가 필요함')
    if qc.get('action')!='created_haplotagged_bam':return {'status':'not_applicable','files':[]}
    bam=Path(qc['haplotagged_bam']);tags=Path(qc['taglist']);index=Path(str(bam)+'.bai')
    for path in [bam,tags,index]:
        if path.is_symlink() or path.parent.resolve()!=folder.resolve() or not path.is_file() or not path.stat().st_size:
            raise ValueError('완성 BAM/태그/인덱스 보존 조건 불일치')
    if qc.get('downstream_alignment')!=str(bam):raise ValueError('후속 입력 BAM 불일치')
    if not qc.get('after',{}).get('hp_present') or not qc.get('taglist_counts',{}).get('total'):
        raise ValueError('HP 태그/태그 목록 검증 근거 없음')
    if bam.stat().st_mtime_ns>qc_json.stat().st_mtime_ns:raise ValueError('QC 이후 완성 BAM이 변경됨')
    subprocess.run([samtools,'quickcheck',str(bam)],check=True,capture_output=True)
    with pysam.AlignmentFile(str(bam),'rb',require_index=True) as source:
        if source.header.to_dict().get('HD',{}).get('SO')!='coordinate' or not source.references:
            raise ValueError('완성 BAM의 좌표 정렬/reference 검증 실패')
        source.get_index_statistics()
    reports=[];parent=folder/'parallel_haplotag'
    if parent.is_symlink():raise ValueError('조각 작업 폴더 symlink 거부')
    for path in parent.glob('*/execution.json'):
        if path.parent.is_symlink() or path.is_symlink():continue
        row=json.loads(path.read_text())
        if row.get('sample')==qc['sample'] and row.get('input_alignment')==qc['input_alignment']:
            reports.append((path,row))
    successful=[(p,r) for p,r in reports if r.get('status')=='success' and r.get('output_bam')==str(bam)
                and float(r.get('finished',0))*1e9<=qc_json.stat().st_mtime_ns]
    if not successful:return {'status':'no_verified_parallel_run','files':[]}
    latest=max(float(r['finished']) for _,r in successful)
    files=[]
    for path,row in reports:
        if row.get('status') not in {'success','failed'} or float(row.get('finished',latest+1))>latest:continue
        for candidate in path.parent.iterdir():
            if not (re.fullmatch(r'\d{4}\.bam',candidate.name) or candidate.name in {'unmapped.bam','combined.bam'}):continue
            if candidate.is_symlink() or not candidate.is_file():continue
            value=candidate.stat()
            if value.st_nlink!=1:continue
            files.append({'path':str(candidate),'fingerprint':fingerprint(candidate),'bytes':value.st_size,'allocated_bytes':value.st_blocks*512})
    result={'schema':'mfm_haplotag_shard_cleanup_v1','sample':qc['sample'],'status':'planned',
            'published_bam':str(bam),'qc_sha256':hashlib.sha256(qc_json.read_bytes()).hexdigest(),
            'files':files,'allocated_bytes':sum(r['allocated_bytes'] for r in files),'time':time.time()}
    if dry_run:return result
    if not files:return dict(result,status='already_clean')
    opened=opened_files([Path(row['path']) for row in files])
    if opened:return dict(result,status='deferred_open_files',opened=opened)
    protected={str(p):fingerprint(p) for p in [bam,tags,index,qc_json]}
    for row in files:
        path=Path(row['path'])
        if path.is_symlink() or fingerprint(path)!=row['fingerprint']:raise ValueError('정리 계획 후 조각 변경됨')
    for row in files:Path(row['path']).unlink()
    if any(fingerprint(Path(path))!=value for path,value in protected.items()):raise ValueError('보존 대상 변경됨')
    result['status']='complete';atomic(parent/'shard_cleanup.json',result)
    return result


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qc-json',required=True,type=Path)
    parser.add_argument('--samtools',default='samtools')
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    print(json.dumps(cleanup_parallel_shards(args.qc_json,args.samtools,dry_run=args.dry_run),ensure_ascii=False))
