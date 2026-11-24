"""NVIDIA Parabricks 호출·장치 공개·host RSS 보호를 공용 모듈로 관리함."""
from pathlib import Path
import json,os,re,resource,shutil,signal,subprocess,tempfile,time,threading

DEFAULT_IMAGE='nvcr.io/nvidia/clara/clara-parabricks:4.8.0-1'

def minimap_gpu_diagnostics(log):
    """GPU 라벨과 실제 seed/chain recovery를 구분해 원시 수치를 기록함."""
    if not log or not Path(log).is_file():return {}
    result={}
    with Path(log).open(errors='replace') as source:
        for line in source:
            match=re.search(r'GPU (seed|chain) kernel time: ([0-9.eE+-]+) s',line)
            if match:result['gpu_'+match[1]+'_kernel_seconds']=float(match[2])
            match=re.search(r'GPU seeding/chaining recovery mode triggered: (\d+)',line)
            if match:result['gpu_seeding_chaining_recovery_count']=int(match[1])
    return result

class Runtime:
    def __init__(self):
        self.kind=os.environ.get('MFM_GPU_RUNTIME','auto')
        self.image=os.environ.get('MFM_PARABRICKS_IMAGE',DEFAULT_IMAGE)
        self.binary=os.environ.get('MFM_PBRUN','pbrun')
        if self.kind=='auto':
            self.kind=next((k for k,tool in [('native',self.binary),('apptainer','apptainer'),
                ('docker','docker')] if shutil.which(tool)),None)

    def validate(self):
        if self.kind not in {'native','docker','apptainer'}:
            raise RuntimeError('Parabricks 실행환경 없음. pbrun 또는 Docker/Apptainer와 GPU 이미지가 필요함')
        tool=self.binary if self.kind=='native' else self.kind
        if not shutil.which(tool):raise RuntimeError('GPU 실행 도구 없음: '+tool)
        if self.kind!='native' and not self.image:raise ValueError('Parabricks 이미지가 필요함')

    def command(self,args,gpu_uuid,paths=(),host_limit=None,cpu_limit=None):
        self.validate()
        if self.kind=='native':return [self.binary,*map(str,args)]
        mapping={str(Path(p).absolute()):str(Path(p).resolve()) for p in paths}
        args=[mapping.get(str(value),str(value)) for value in args]
        binds=set()
        for item in paths:
            path=Path(item).expanduser().absolute()
            anchor=path if path.is_dir() else path.parent
            binds.add(str(anchor.resolve()))
        binds=sorted(p for p in binds if not any(p!=x and Path(x) in Path(p).parents for x in binds))
        if self.kind=='docker':
            command=['docker','run','--rm','--user',f'{os.getuid()}:{os.getgid()}',
                     '--env','CUDA_VISIBLE_DEVICES='+(','.join(str(i) for i in range(len(gpu_uuid.split(',')))) if gpu_uuid else ''),'--env','NVIDIA_VISIBLE_DEVICES='+(gpu_uuid or 'none'),
                     '--env','HOME=/tmp','--env','TMPDIR=/tmp',
                     '--memory',str(int(host_limit or os.environ.get('MFM_GPU_HOST_MEMORY_BYTES',str(128*1024**3)))),
                     '--cpus',str(max(1,int(cpu_limit or os.environ.get('MFM_GPU_HOST_CPU_LIMIT','16'))))]
            if gpu_uuid:command+=['--gpus',('\"device='+gpu_uuid+'\"') if ',' in gpu_uuid else 'device='+gpu_uuid]
            for path in binds:command+=['--volume',path+':'+path]
            return [*command,self.image,'pbrun',*map(str,args)]
        image=self.image if Path(self.image).is_file() or self.image.startswith(('docker://','oras://')) else 'docker://'+self.image
        command=['apptainer','exec','--nv']
        for path in binds:command+=['--bind',path+':'+path]
        return [*command,image,'pbrun',*map(str,args)]

    def run(self,args,gpu,paths=(),host_limit=None,log=None,cpu_limit=None,check=None):
        from .memory_policy import process_rss
        host_limit=int(host_limit or os.environ.get('MFM_GPU_HOST_MEMORY_BYTES',str(128*1024**3)))
        gpu_uuid=gpu.device['uuid'] if gpu is not None else None
        command=self.command(args,gpu_uuid,paths,host_limit=host_limit,cpu_limit=cpu_limit)
        env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES=gpu_uuid or '',NVIDIA_VISIBLE_DEVICES=gpu_uuid or 'none')
        if self.kind=='apptainer':env['APPTAINERENV_CUDA_VISIBLE_DEVICES']=gpu_uuid or ''
        soft,hard=resource.getrlimit(resource.RLIMIT_AS)
        def cuda_child():
            # CUDA의 큰 가상 주소 예약은 RSS와 다름. 아래 실제 RSS 검사로 host 상한을 보호함.
            if hard==resource.RLIM_INFINITY:resource.setrlimit(resource.RLIMIT_AS,(resource.RLIM_INFINITY,hard))
        started=time.time();peak=0
        stream=open(log,'ab') if log else None
        handlers={}
        def interrupted(signum,frame):
            raise SystemExit(128+signum)
        try:
            process=subprocess.Popen(command,env=env,stdout=stream or subprocess.DEVNULL,
                stderr=stream,preexec_fn=cuda_child,start_new_session=True)
            try:
                if threading.current_thread() is threading.main_thread():
                    for sig in [signal.SIGTERM,signal.SIGINT]:
                        handlers[sig]=signal.signal(sig,interrupted)
                while process.poll() is None:
                    if check is not None:check()
                    rss=process_rss({'gpu':process.pid})['gpu'];peak=max(peak,rss)
                    if rss>host_limit:
                        os.killpg(process.pid,signal.SIGTERM);process.wait()
                        raise MemoryError(f'GPU 도구 host RSS 예약 초과: {rss}/{host_limit}')
                    time.sleep(.5)
                if process.returncode:raise subprocess.CalledProcessError(process.returncode,command)
            except BaseException:
                if process.poll() is None:os.killpg(process.pid,signal.SIGTERM);process.wait()
                raise
        finally:
            for sig,handler in handlers.items():signal.signal(sig,handler)
            if stream:stream.close()
        return {**minimap_gpu_diagnostics(log),'runtime':self.kind,'image':self.image if self.kind!='native' else None,
                'command':command,'gpu_uuid':gpu_uuid,'gpu_uuids':gpu_uuid.split(',') if gpu_uuid else [],'gpu_count':len(gpu_uuid.split(',')) if gpu_uuid else 0,'gpu_reserved_bytes':gpu.reserved if gpu is not None else 0,
                'started_at':started,'elapsed_seconds':time.time()-started,'peak_host_rss_bytes':(None if self.kind=='docker' else peak),
                'host_rss_limit_bytes':host_limit,'rss_scope':('Docker client process tree; workload is protected by container memory limit' if self.kind=='docker' else 'GPU tool process tree'),
                'container_memory_limit_bytes':host_limit if self.kind=='docker' else None,'address_space_policy':'CUDA virtual address space; per-process-tree RSS guard'}
