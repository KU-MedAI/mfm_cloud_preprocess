"""선택된 샘플만 다운로드하고 크기/checksum/소유권을 기록한다."""
from __future__ import annotations
import fcntl, hashlib, json, os, shutil, subprocess, time
from functools import lru_cache
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request,urlopen
from .planner import read_tsv
from .download_telemetry import TransferLog

def safe_path(root, relative):
    root=Path(root).resolve();rel=Path(relative)
    if rel.is_absolute() or '..' in rel.parts:raise ValueError('상대 경로 밖으로 나갈 수 없습니다: '+relative)
    p=root/rel
    if root not in p.resolve().parents:raise ValueError('입력 저장소 밖 경로: '+str(p))
    return p

def verify(path,size=None,checksum=None,algorithm='md5'):
    p=Path(path)
    if not p.is_file() or p.stat().st_size==0:return False
    if size and p.stat().st_size!=int(size):return False
    if checksum:
        digest=hashlib.new(algorithm)
        with p.open('rb') as f:
            while b:=f.read(8<<20):digest.update(b)
        if digest.hexdigest().lower()!=checksum.lower():return False
    return True

def fetch(resource, data_root, state_root):
    target=safe_path(data_root,resource['path']);target.parent.mkdir(parents=True,exist_ok=True)
    # 입력 저장소별 잠금: 서로 다른 배치 state_root라도 같은 파일을 동시에 쓰지 않는다.
    state=Path(data_root)/'.download_records';state.mkdir(parents=True,exist_ok=True)
    key=hashlib.sha256(str(target).encode()).hexdigest();record=state/(key+'.json')
    size=resource.get('size_bytes') or None
    checksum=resource.get('checksum') or None;algorithm=resource.get('checksum_type','md5')
    with (state/(key+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        prior=json.loads(record.read_text()) if record.exists() else {}
        signature={'size':target.stat().st_size,'mtime_ns':target.stat().st_mtime_ns,'ctime_ns':target.stat().st_ctime_ns} if target.is_file() else None
        generated=prior.get('status')=='generated_index'
        cached=signature is not None and prior.get('status') in {'verified','generated_index'} and prior.get('url')==resource['url'] and prior.get('signature')==signature
        if generated:
            cached=cached and bool(prior.get('generated_from'))
            if cached:
                try:
                    s=Path(prior['generated_from']).stat()
                    cached=prior.get('generated_from_signature')=={'size':s.st_size,'mtime_ns':s.st_mtime_ns,'ctime_ns':s.st_ctime_ns}
                except OSError:
                    cached=False
        else:
            cached=cached and (not size or signature['size']==int(size)) and (prior.get('checksum') or '').lower()==(checksum or '').lower() and (not checksum or prior.get('checksum_type','md5')==algorithm)
        if cached:
            return target
        if not generated and verify(target,size,checksum,algorithm):
            s=target.stat()
            after={'size':s.st_size,'mtime_ns':s.st_mtime_ns,'ctime_ns':s.st_ctime_ns}
            if signature!=after:raise ValueError('검증 중 입력 파일이 변경되었습니다: '+str(target))
            prior.update(path=str(target),url=resource['url'],owned=bool(prior.get('owned',False)),status='verified',
                         size_bytes=s.st_size,checksum=checksum,checksum_type=algorithm,signature=after)
            record.write_text(json.dumps(prior,ensure_ascii=False)+'\n')
            return target
        if target.exists() and not prior.get('owned'):raise ValueError('기존 파일 검증 실패: 자동으로 덮어쓰지 않습니다: '+str(target))
        partial=Path(str(target)+'.part')
        # curl은 재시도와 HTTP Range 재개를 지원한다. 인증 정보는 전달본에 포함하지 않는다.
        url=resource['url'];scheme=urlparse(url).scheme
        with TransferLog(url,target,partial,state_root,data_root,kind=resource.get('kind','')) as telemetry:
            with telemetry.phase('transfer'):
                if scheme=='file':
                    shutil.copyfile(Path(urlparse(url).path),partial)
                elif scheme in {'https','http','ftp'}:
                    from .download_retry import transfer
                    transfer(url,partial,lambda p:bool(size or checksum) and verify(p,size,checksum,algorithm))
                else:raise ValueError('지원하지 않는 다운로드 URL: '+url)
                telemetry.progress(partial.stat().st_size)
            with telemetry.phase('verification'):
                if not verify(partial,size,checksum,algorithm):raise ValueError('다운로드 크기/checksum 검증 실패: '+str(target))
                checked=None
                if resource.get('kind') in {'wgs_alignment','wgs_read','rna_alignment'} and target.suffix in {'.bam','.cram'}:
                    from .io_integrity import validate
                    checked=validate(partial,state/(key+'.integrity.json'))
                os.replace(partial,target)
                if checked is not None:
                    from .integrity_cache import publish_relocated
                    publish_relocated(checked,target,state/(key+'.integrity.json'))
        record.write_text(json.dumps({'path':str(target),'url':url,'owned':True,'status':'verified',
          'size_bytes':target.stat().st_size,'checksum':checksum,'checksum_type':algorithm,'time':time.time(),
          'transfer_metrics':telemetry.summary(),
          'signature':{'size':target.stat().st_size,'mtime_ns':target.stat().st_mtime_ns,'ctime_ns':target.stat().st_ctime_ns}},ensure_ascii=False)+'\n')
    return target

def fetch_or_index(resource,data_root,state_root,threads=1):
    try:
        return fetch(resource,data_root,state_root)
    except subprocess.CalledProcessError as exc:
        if resource.get('kind') not in {'wgs_index','rna_index','sample_tbi','shared_tbi','illumina_panel_index'}:raise
        # 원본 다운로드 성공 후에도 공개 index URL이403/404일 수 있다. alignment/VCF로 직접 재생성한다.
        target=safe_path(data_root,resource['path'])
        base=Path(str(target).removesuffix('.bai').removesuffix('.crai').removesuffix('.tbi'))
        if not base.is_file():raise ValueError('index 생성용 원본이 없습니다: '+str(base)) from exc
        state=Path(data_root)/'.download_records';state.mkdir(parents=True,exist_ok=True);key=hashlib.sha256(str(target).encode()).hexdigest()
        record=state/(key+'.json')
        with (state/(key+'.lock')).open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            old=json.loads(record.read_text()) if record.exists() else {}
            if target.exists() and not old.get('owned'):raise ValueError('기존 사용자 index는 덮어쓰지 않습니다: '+str(target))
            if target.suffix=='.tbi':command=['bcftools','index','--tbi','--force','--output',str(target),str(base)]
            else:command=['samtools','index','-@',str(max(1,threads)),'-o',str(target),str(base)]
            before=base.stat();subprocess.run(command,check=True);after=base.stat()
            if not target.is_file() or not target.stat().st_size or (before.st_size,before.st_mtime_ns,before.st_ctime_ns)!=(after.st_size,after.st_mtime_ns,after.st_ctime_ns):
                raise ValueError('index 생성 중 원본이 변경되었거나 출력이 없습니다.')
            st=target.stat()
            record.write_text(json.dumps({'path':str(target),'url':resource['url'],'owned':True,'status':'generated_index',
              'command':command,'reason':'원본 index 다운로드 실패: '+str(exc),
              'generated_from':str(base),'generated_from_signature':{'size':after.st_size,'mtime_ns':after.st_mtime_ns,'ctime_ns':after.st_ctime_ns},
              'signature':{'size':st.st_size,'mtime_ns':st.st_mtime_ns,'ctime_ns':st.st_ctime_ns}},ensure_ascii=False,indent=2)+'\n')
            resource['index_fallback']='원본으로 로컬 index 생성';resource['actual_index_bytes']=st.st_size
        return target

@lru_cache(maxsize=8)
def resource_tables(root):
    root=Path(root);cap=root/'manifests/capacity_catalog_260929_complete'
    records={r['resource_id']:r for r in read_tsv(cap/'raw_resource_inventory.tsv')}
    links={}
    for r in read_tsv(cap/'raw_resource_links.tsv'):
        links.setdefault((r['consumer_type'],r['consumer_id']),set()).add(r['resource_id'])
    return records,links

def sample_resources(root,row):
    records,links=resource_tables(str(Path(root).resolve()))
    ids=links.get(('sample',row['sample_id']),set())|links.get(('project',row['project_id']),set())
    out=[]
    for key in sorted(ids):
        r=records[key];kind=r['resource_kind'];rel=r['relative_path']
        if kind.startswith('reference'):continue # 참조 준비 모듈에서 검증된 역할별 URL 사용
        if kind in {'sample_tsv'}:continue
        if kind in {'wgs_index','rna_index'}:
            # 기존 용량 목록의 index relative_path가 alignment 경로와 같은 행을 보정한다.
            suffix='.crai' if '.cram' in urlparse(r['source_url']).path else '.bai'
            if not rel.endswith(suffix):rel+=suffix
        if not rel:
            base=Path(urlparse(r['source_url']).path).name
            if kind.startswith('illumina_panel'):rel='common_vcf/'+base
            elif kind.startswith('shared') or kind.startswith('sample_'):rel='public_vcf/'+row['project_id']+'/'+base
            else:raise ValueError('다운로드 경로 누락: '+kind)
        algorithm=r['checksum_type'];checksum=r['checksum_value']
        if algorithm not in {'md5','sha256'}:algorithm='md5';checksum=''
        out.append({'path':rel,'url':r['source_url'],'size_bytes':r['size_bytes'],
           'checksum_type':algorithm,'checksum':checksum,'shared':r['scope']=='shared','kind':kind})
    return out

def remove_owned(resources,data_root,state_root):
    removed=[]
    for r in resources:
        if r.get('shared'):continue
        # 이 함수는 read 원본만 정리함. 비공유 개인 VCF·외부 작업 VCF는 completion_cleanup에서 완료 후 정리함.
        if r.get('kind') not in {'wgs_alignment','wgs_read','wgs_index','rna_alignment','rna_index'}:continue
        p=safe_path(data_root,r['path']);key=hashlib.sha256(str(p).encode()).hexdigest()
        state=Path(data_root)/'.download_records';record=state/(key+'.json')
        if not record.exists():continue
        with (state/(key+'.lock')).open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX);d=json.loads(record.read_text())
            if d.get('owned') and d.get('path')==str(p):
                if p.exists():p.unlink();removed.append(str(p))
                d['status']='removed_after_package_pass';record.write_text(json.dumps(d,ensure_ascii=False)+'\n')
    return removed
