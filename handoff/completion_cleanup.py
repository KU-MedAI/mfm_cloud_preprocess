"""최종검증후샘플외부의전용원본·파생VCF만제거함."""
from pathlib import Path
import os,re,json,hashlib,fcntl,time
from tools.compact_manifest import completed_contigs

def guard(root,sample,mode):
    from .cleanup import cleanup_complete
    root=Path(root).resolve()
    if root.name!=sample or not cleanup_complete(root,sample,mode,external=False):raise ValueError('정리후PASS 완료샘플만외부정리가능함')
    report=json.loads((root/'qc/final_package_validation.json').read_text())
    if report.get('sample')!=sample or report.get('mode')!=mode or report.get('status')!='pass' or report.get('failed_checks'):raise ValueError('최종QC불일치')
    delivery=json.loads((root/'qc/delivery_final_validation.json').read_text())
    if delivery.get('status')!='pass' or delivery.get('errors'):raise ValueError('정리후실제검증보고서불일치')
    manifests=list((root/'compact').glob('*/manifest.json'))
    if len(manifests)!=1:raise ValueError('compact manifest개수불일치')
    m=json.loads(manifests[0].read_text())
    if m.get('status')!='complete' or m.get('validation_level')!='full':raise ValueError('full compact필요')
    expected=completed_contigs(m,manifests[0].parent,sample)
    traces=list((root/'qc/stage12_trace').glob('*/trace_manifest.json'))
    if not expected or len(traces)!=1:raise ValueError('독립trace누락')
    t=json.loads(traces[0].read_text())
    if t.get('status')!='pass' or {c for c,v in t.get('contigs',{}).items() if v=='pass'}!=expected or Path(t.get('compact_root','')).resolve()!=manifests[0].parent.resolve():raise ValueError('독립trace불일치')
    for c in expected:
        required=['hap_sequence.zarr','readset.HAP1.zarr','readset.HAP2.zarr','depth.zarr','windows.tsv']
        if m.get('wgs_payload')=='sequence_only':
            if report.get('wgs_payload')!='sequence_only':raise ValueError('검증된 WGS payload 불일치')
            required=['hap_sequence.zarr','windows.tsv']
        if m.get('wgs_payload')=='pileup_feature':
            if report.get('wgs_payload')!='pileup_feature':raise ValueError('검증된 pileup payload 불일치')
            required=['hap_sequence.zarr','pileup.zarr','windows.tsv']
        for name in required:
            if not (manifests[0].parent/c/name).exists():raise ValueError('필수모델입력누락')
    return {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'qc/final_package_validation.json',root/'qc/delivery_final_validation.json',root/'metadata/cleanup_completed.json',manifests[0]]}

def _fingerprint(p):
    s=p.stat();return [s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns]

def plan(root,sample,mode,data_root,output_root,allow_unowned=False):
    root=Path(root).resolve();data=Path(data_root).resolve();output=Path(output_root).resolve()
    if output not in root.parents or root==data or root in data.parents:raise ValueError('외부정리경로불일치')
    marker=data/'.mfm_handoff_store.json'
    if not marker.is_file() or json.loads(marker.read_text()).get('schema')!='mfm_handoff_store_v1':raise ValueError('관리원본저장소표식필요')
    hashes=guard(root,sample,mode);meta=json.loads((root/'metadata/sample.json').read_text());source=meta.get('source_sample_id') or sample
    inputs=root/'metadata/download_inputs.json'
    shared={str((data/r['path']).resolve()) for r in json.loads(inputs.read_text()) if r.get('shared')} if inputs.exists() else set()
    prefix=re.compile(r'^(?:'+re.escape(sample)+'|'+re.escape(source)+r')(?=[._-]|$)');candidates={}
    def add(p,kind):
        if p.is_symlink():raise ValueError('외부정리symlink거부: '+str(p))
        if not p.is_file():return
        real=p.resolve()
        if root==real or root in real.parents or not (data in real.parents or output/'_work' in real.parents):raise ValueError('허용범위밖삭제거부')
        if real in candidates or str(real) in shared:return
        if kind in {'raw','private_vcf'} and not allow_unowned:
            key=hashlib.sha256(str(real).encode()).hexdigest();record=data/'.download_records'/(key+'.json')
            if not record.is_file() or not json.loads(record.read_text()).get('owned'):return
        candidates[real]={'path':str(real),'kind':kind,'fingerprint':_fingerprint(real),'bytes':real.stat().st_size}
    # 이번 실행에 사용하지 않은 기존 WGS/RNA 원본은 후속 처리를 위해 보존함.
    raw_areas=[] if meta.get('rna_processing_status')=='deferred' else ['wgs','rna']
    for area in raw_areas:
        for p in (data/area).glob('*/*'):
            if p.name not in {sample,source}:continue
            if p.is_symlink():raise ValueError('원본샘플디렉터리symlink거부')
            for f in p.rglob('*'):
                if f.is_dir():continue
                if f.name.endswith(('.json','.lock')) and not f.name.endswith('.part.http.json'):continue
                add(f,'raw')
    # 샘플명이없는공용VCF·reference는대상이아님.
    for f in (data/'public_vcf').glob('*/*'):
        if prefix.match(f.name) and str(f.resolve()) not in shared:add(f,'private_vcf')
    work=output/'_work'
    if work.exists():
        for here,dirs,files in os.walk(work,followlinks=False):
            if 'persample_vcf' not in Path(here).parts:continue
            for n in files:
                if prefix.match(n) and n.endswith(('.vcf','.vcf.gz','.tbi','.csi','.vcf.gz.part','.vcf.gz.tmp')):add(Path(here)/n,'work_vcf')
    return {'schema':'mfm_external_cleanup_v1','sample':sample,'mode':mode,'root':str(root),'data_root':str(data),'output_root':str(output),'model_guard_hashes':hashes,'files':list(candidates.values()),'logical_bytes':sum(r['bytes'] for r in candidates.values()),'allow_unowned':allow_unowned}

def _opened(paths):
    found=[]
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:
            for fd in (p/'fd').iterdir():
                try:target=os.readlink(fd)
                except OSError:continue
                if target in paths or (target.startswith('/') and str(Path(target).resolve()) in paths):found.append({'pid':int(p.name),'path':target})
        except (PermissionError,FileNotFoundError,ProcessLookupError):continue
    return found

def cleanup_external(root,sample,mode,data_root,output_root,allow_unowned=False,dry_run=False):
    root=Path(root).resolve()
    with (root/'metadata/.external_cleanup.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        probe=root/'metadata/external_cleanup_in_progress.json'
        if not dry_run:
            tmp=probe.with_suffix('.json.tmp');tmp.write_text(json.dumps({'schema':'mfm_external_cleanup_v1','sample':sample,'mode':mode,'status':'checking','time':time.time()})+'\n');os.replace(tmp,probe)
        try:
            result=plan(root,sample,mode,data_root,output_root,allow_unowned)
        except Exception as error:
            if not dry_run:
                tmp=probe.with_suffix('.json.tmp');tmp.write_text(json.dumps({'schema':'mfm_external_cleanup_v1','sample':sample,'mode':mode,'status':'failed','error':str(error),'time':time.time()})+'\n');os.replace(tmp,probe)
            raise
        result['time']=time.time();result['status']='planned';result['removed']=[]
        if dry_run:return result
        def save(name):
            p=root/'metadata'/name;tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');os.replace(tmp,p)
        result['status']='running';save('external_cleanup_in_progress.json')
        try:
            opened=_opened({f['path'] for f in result['files']})
            if opened:raise ValueError('삭제대상열린파일존재: '+str(opened[:4]))
            for f in result['files']:
                p=Path(f['path'])
                if p.is_symlink() or _fingerprint(p)!=f['fingerprint']:raise ValueError('계획후변경파일삭제거부: '+str(p))
                p.unlink();result['removed'].append(str(p))
            if any(Path(f['path']).exists() for f in result['files']):raise ValueError('외부중간잔존')
            if any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in result['model_guard_hashes'].items()):raise ValueError('모델검증해시변경')
            result['status']='complete';save('external_cleanup_completed.json');(root/'metadata/external_cleanup_in_progress.json').unlink(missing_ok=True)
        except Exception as e:
            result['status']='failed';result['error']=str(e);save('external_cleanup_in_progress.json');raise
        return result
