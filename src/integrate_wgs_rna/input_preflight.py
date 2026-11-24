"""공통 reference M5와 VCF REF 전수 검사를 병렬화하고 동일 입력 성공만 재사용함."""
import argparse,fcntl,hashlib,json,os,subprocess,tempfile,time
from pathlib import Path
import pysam
from .input_parallel import cpu_budget,parallel_calls
from .memory_policy import GIB
SCHEMA='mfm_input_preflight_v1'


def identity(path):
    p=Path(path).resolve();s=p.stat()
    if not p.is_file() or not s.st_size:raise ValueError('검사 파일 없음/비어 있음: '+str(p))
    return {'path':str(p),'bytes':s.st_size,'mtime_ns':s.st_mtime_ns,'ctime_ns':s.st_ctime_ns,'inode':s.st_ino}


def key(value):return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()


def atomic(path,content):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False) as f:f.write(content);tmp=Path(f.name)
    try:os.replace(tmp,path)
    finally:tmp.unlink(missing_ok=True)


def _m5(reference,contig):
    with pysam.FastaFile(str(reference)) as f:
        length=f.get_reference_length(contig);digest=hashlib.md5()
        for start in range(0,length,8*1024**2):digest.update(f.fetch(contig,start,min(length,start+8*1024**2)).encode('ascii'))
    return contig,length,digest.hexdigest()


def reference_m5(reference,contigs,output,cache_root,jobs):
    reference=Path(reference);before={'schema':SCHEMA,'kind':'m5_case_preserved','reference':identity(reference),'index':identity(str(reference)+'.fai'),'contigs':contigs}
    cache=Path(cache_root)/('m5_'+key(before)+'.json');cache.parent.mkdir(parents=True,exist_ok=True)
    with cache.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);cached=False
        if cache.exists():
            try:
                saved=json.loads(cache.read_text());rows=saved['rows'];cached=saved.get('identity')==before and len(rows)==len(contigs) and [r[0] for r in rows]==contigs
            except (ValueError,KeyError,TypeError):pass
        if not cached:
            rows=parallel_calls(_m5,[(reference,c) for c in contigs],min(jobs,len(contigs)),memory_per_task=256*1024**2)
            if before['reference']!=identity(reference) or before['index']!=identity(str(reference)+'.fai'):raise RuntimeError('M5 검사 중 reference 변경됨')
            atomic(cache,json.dumps({'identity':before,'rows':rows},ensure_ascii=False)+'\n')
        atomic(output,''.join(f'{c}\t{length}\t{digest}\n' for c,length,digest in rows))
    return cached


def _validate_vcf(vcf,label,reference,regions,bcftools,cache_root,error_dir,sample):
    started=time.monotonic();cache=None
    err=Path(error_dir)/(sample+'.'+label+'.vcf_ref.err')
    try:
        index=next((Path(str(vcf)+ext) for ext in ['.tbi','.csi'] if Path(str(vcf)+ext).is_file()),None)
        if index is None:raise ValueError('VCF index 없음: '+str(vcf))
        before={'schema':SCHEMA,'kind':'vcf_norm_check_ref_e','vcf':identity(vcf),'index':identity(index),
                'reference':identity(reference),'reference_index':identity(str(reference)+'.fai'),'regions':regions,
                'bcftools':identity(bcftools),'options':['view','-G','-r',regions,'-Ou','norm','--check-ref','e','--fasta-ref','-Ou']}
        cache=Path(cache_root)/('vcf_'+key(before)+'.json');cache.parent.mkdir(parents=True,exist_ok=True)
        with cache.with_suffix('.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            if cache.exists():
                try:
                    saved=json.loads(cache.read_text())
                    if saved.get('status')=='pass' and saved.get('identity')==before:
                        return label,True,f'{label} shared_cache={cache}',True,time.monotonic()-started
                except (OSError,ValueError):pass
            err.parent.mkdir(parents=True,exist_ok=True)
            with err.open('wb') as log:
                view=subprocess.Popen([str(bcftools),'view','-G','-r',regions,'-Ou',str(vcf)],stdout=subprocess.PIPE,stderr=log)
                try:
                    assert view.stdout is not None
                    norm=subprocess.Popen([str(bcftools),'norm','--check-ref','e','--fasta-ref',str(reference),'-Ou','-'],stdin=view.stdout,stdout=subprocess.DEVNULL,stderr=log)
                except BaseException:
                    view.terminate();view.wait();raise
                finally:
                    if view.stdout:view.stdout.close()
                try:rc=norm.wait();view_rc=view.wait()
                except BaseException:
                    for child in [norm,view]:
                        if child.poll() is None:child.terminate()
                    norm.wait();view.wait();raise
            if rc or view_rc:raise ValueError(err.read_text(errors='replace')[-1000:] or f'bcftools 실패 {view_rc}/{rc}')
            for name,path in [('vcf',vcf),('index',index),('reference',reference),('reference_index',str(reference)+'.fai'),('bcftools',bcftools)]:
                if before[name]!=identity(path):raise RuntimeError('REF 검사 중 입력/도구 변경됨: '+str(path))
            atomic(cache,json.dumps({'status':'pass','identity':before,'validated_at':time.time()},ensure_ascii=False)+'\n')
            return label,True,f'{label} autosomes={regions} shared_cache={cache}',False,time.monotonic()-started
    except Exception as error:
        return label,False,str(error).replace('\n',' ').replace('\t',' ')[:1000],False,time.monotonic()-started


def validate_vcfs(items,reference,regions,bcftools,cache_root,error_dir,sample,jobs):
    return parallel_calls(_validate_vcf,[(vcf,label,reference,regions,bcftools,cache_root,error_dir,sample) for vcf,label in items],min(jobs,len(items)),cpu_cost=min(2,cpu_budget()),memory_per_task=2*GIB)


def main():
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='mode',required=True)
    m=sub.add_parser('m5');m.add_argument('--contigs',required=True)
    v=sub.add_parser('vcf');v.add_argument('--items',required=True);v.add_argument('--regions',required=True);v.add_argument('--bcftools',required=True);v.add_argument('--error-dir',required=True);v.add_argument('--sample',required=True)
    for cmd in [m,v]:
        cmd.add_argument('--reference',required=True);cmd.add_argument('--cache-root',required=True);cmd.add_argument('--output',required=True);cmd.add_argument('--jobs',type=int,default=cpu_budget())
    a=p.parse_args()
    if a.jobs<1:raise ValueError('jobs는 양수여야 함')
    if a.mode=='m5':reference_m5(a.reference,a.contigs.split(','),a.output,a.cache_root,a.jobs)
    else:
        items=[line.split('\t') for line in Path(a.items).read_text().splitlines() if line]
        rows=validate_vcfs(items,a.reference,a.regions,a.bcftools,a.cache_root,a.error_dir,a.sample,min(a.jobs,max(1,cpu_budget()//2)))
        atomic(a.output,''.join(f'{label}\t{int(ok)}\t{detail}\t{int(cached)}\t{seconds:.6f}\n' for label,ok,detail,cached,seconds in rows))
        print('[PREFLIGHT-VCF] '+json.dumps({'tasks':len(rows),'pass':sum(x[1] for x in rows),'cached':sum(x[3] for x in rows),'jobs':min(a.jobs,max(1,cpu_budget()//2))}),flush=True)


if __name__=='__main__':main()
