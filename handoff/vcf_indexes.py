"""검증전 누락된 public VCF index만 잠금·원자교체로생성함."""
from pathlib import Path
import argparse,fcntl,json,os,subprocess,tempfile,time

def ensure_index(source,bcftools='bcftools',threads=1):
    source=Path(source)
    if not source.is_file():raise FileNotFoundError(source)
    with Path(str(source)+'.index.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if any(Path(str(source)+s).is_file() and Path(str(source)+s).stat().st_size>0 for s in ['.tbi','.csi']):return 'existing'
        target=Path(str(source)+'.csi');before=source.stat()
        fd,name=tempfile.mkstemp(prefix='.'+target.name+'.',dir=target.parent);os.close(fd);temporary=Path(name)
        try:
            subprocess.run([bcftools,'index','--csi','--threads',str(max(0,int(threads)-1)),'--force','--output',str(temporary),str(source)],check=True)
            after=source.stat()
            if (before.st_size,before.st_mtime_ns,before.st_ctime_ns)!=(after.st_size,after.st_mtime_ns,after.st_ctime_ns):raise ValueError('index 생성중 원본VCF가변경됨')
            import pysam
            with pysam.VariantFile(str(source),index_filename=str(temporary)) as reader:
                next(reader.fetch(next(iter(reader.header.contigs))),None)
            if not temporary.stat().st_size:raise ValueError('생성index가비어있음')
            os.replace(temporary,target)
            print('[PUBLIC-VCF-INDEX] '+json.dumps({'source':str(source),'index':str(target),'threads':threads,'created_at':time.time()},ensure_ascii=False),flush=True)
            return 'created'
        finally:temporary.unlink(missing_ok=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--bcftools',default='bcftools');parser.add_argument('sources',nargs='+');args=parser.parse_args()
    caps=[4]
    for key in ['MFM_CPU_POOL_GRANTED','MFM_SAMPLE_CPU_BUDGET','SLURM_CPUS_PER_TASK']:
        if os.environ.get(key):caps.append(int(os.environ[key]))
    if os.environ.get('MFM_CPU_POOL_LEASE_FILE'):caps.append(int(json.loads(Path(os.environ['MFM_CPU_POOL_LEASE_FILE']).read_text())['cpus']))
    for source in args.sources:ensure_index(source,args.bcftools,max(1,min(caps)))

if __name__=='__main__':main()
