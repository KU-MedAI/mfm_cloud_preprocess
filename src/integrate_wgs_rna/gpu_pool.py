"""허용 GPU 목록과 VRAM 예약을 flock으로 공유함. 다른 학습 GPU는 선택하지 않음."""
from contextlib import contextmanager
from pathlib import Path
import atexit,csv,fcntl,json,os,subprocess,tempfile,time,uuid

GIB=1024**3

def inventory():
    result=subprocess.run(['nvidia-smi','--query-gpu=index,uuid,pci.bus_id,memory.total,memory.used,utilization.gpu',
                           '--format=csv,noheader,nounits'],text=True,capture_output=True,check=True)
    devices=[dict(index=int(v[0]),uuid=v[1].strip(),pci=v[2].strip(),
                 total=int(v[3])*1024**2,used=int(v[4])*1024**2,utilization=int(v[5]))
            for v in csv.reader(result.stdout.splitlines())]
    # 드라이버 잔여 메모리만 있는 장치는 실제 compute PID까지 확인함.
    try:
        apps=subprocess.run(['nvidia-smi','--query-compute-apps=gpu_uuid,pid',
                             '--format=csv,noheader,nounits'],text=True,capture_output=True,check=True)
        pids={d['uuid']:[] for d in devices}
        for row in csv.reader(apps.stdout.splitlines()):
            if len(row)>=2 and row[1].strip().isdigit() and row[0].strip() in pids:pids[row[0].strip()].append(int(row[1]))
        for d in devices:d['compute_pids']=pids[d['uuid']]
    except subprocess.CalledProcessError:pass
    return devices

def allowed_devices(devices, requested=None, environment=None):
    env=os.environ if environment is None else environment
    selected=list(devices)
    # CUDA/컨테이너/Slurm이 공개한 장치 범위를 CLI가 넓힐 수 없음.
    for raw in [env.get('NVIDIA_VISIBLE_DEVICES'),env.get('CUDA_VISIBLE_DEVICES'),
                requested if requested is not None else env.get('MFM_GPU_DEVICES','auto')]:
        if raw is None or str(raw).strip() in {'auto','all'}:continue
        tokens={v.strip() for v in str(raw).split(',') if v.strip()}
        selected=[d for d in selected if str(d['index']) in tokens or any(
                  t.startswith('GPU-') and d['uuid'].startswith(t) for t in tokens)]
    return selected

def _alive(row):
    try:
        text=(Path('/proc')/str(row['pid'])/'stat').read_text()
        fields=text[text.rfind(')')+2:].split()
        return fields[0] not in {'Z','X'} and int(fields[19])==int(row['start_ticks'])
    except FileNotFoundError:return False

def _ticks():
    text=Path('/proc/self/stat').read_text();return int(text[text.rfind(')')+2:].split()[19])

def _save(path,value):
    p=path.with_name('.'+path.name+'.'+str(os.getpid())+'.part')
    p.write_text(json.dumps(value,ensure_ascii=False)+'\n');os.replace(p,path)

class Lease:
    def __init__(self,device,path,reserved):self.device,self.path,self.reserved=device,path,reserved
    def close(self):
        with (self.path.parent.parent/'pool.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX);self.path.unlink(missing_ok=True)

class LeaseGroup:
    def __init__(self, leases):
        self.leases=leases;self.devices=[lease.device for lease in leases]
        self.device={'uuid':','.join(d['uuid'] for d in self.devices)}
        self.reserved=sum(lease.reserved for lease in leases)
    def close(self):
        with (self.leases[0].path.parent.parent/'pool.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            for lease in self.leases:lease.path.unlink(missing_ok=True)


def acquire(required_bytes, *, purpose, exclusive=False, timeout=None, count=1):
    if isinstance(required_bytes,bool) or int(required_bytes)<=0:raise ValueError('GPU VRAM 예약은 양수여야 함')
    required_bytes=int(required_bytes)
    if isinstance(count,bool) or not isinstance(count,int) or count<1:raise ValueError('GPU 개수는 양의 정수여야 함')
    if count>1 and not exclusive:raise ValueError('다중GPU는 독점 원자 예약만 지원함')
    root=Path(os.environ.get('MFM_GPU_POOL_DIR',str(Path(tempfile.gettempdir())/f'mfm_gpu_pool_{os.getuid()}')))
    directory=root/'leases';directory.mkdir(parents=True,exist_ok=True)
    deadline=time.monotonic()+float(os.environ.get('MFM_GPU_WAIT_SECONDS','300') if timeout is None else timeout)
    while True:
        devices=allowed_devices(inventory())
        if len(devices)<count:raise RuntimeError('요청 GPU 개수가 허용된 장치 수보다 큼. --gpu-devices 및 CUDA_VISIBLE_DEVICES 확인 필요함')
        with (root/'pool.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            active=[]
            for path in directory.glob('*.json'):
                row=json.loads(path.read_text())
                if _alive(row):active.append(row)
                else:path.unlink(missing_ok=True)
            options=[]
            for device in devices:
                owned=[r for r in active if r['gpu_uuid']==device['uuid']]
                if any(r['exclusive'] for r in owned) or (exclusive and owned):continue
                # 자체 예약이 없는 GPU에 사용 흔적이 있으면 다른 작업용으로 보존함.
                if not owned and device.get('compute_pids'):continue
                if not owned and (device['used']>0 or device['utilization']>0):
                    idle_driver=(device['used']<=16*(1<<20) and device['utilization']==0 and device.get('compute_pids')==[])
                    if not idle_driver:continue
                reserved=sum(int(r['reserved_bytes']) for r in owned)
                margin=max(GIB,device['total']//20)
                if max(reserved,device['used'])+required_bytes>device['total']-margin:continue
                options.append((reserved/device['total'],len(owned),device['used']/device['total'],device['index'],device))
            if len(options)>=count:
                leases=[]
                try:
                    for option in sorted(options,key=lambda v:v[:-1])[:count]:
                        device=option[-1];path=directory/(uuid.uuid4().hex+'.json')
                        value=dict(schema='mfm_gpu_lease_v1',pid=os.getpid(),start_ticks=_ticks(),
                            gpu_uuid=device['uuid'],gpu_index=device['index'],reserved_bytes=(device['total'] if exclusive else required_bytes),
                            exclusive=bool(exclusive),purpose=purpose,sample=os.environ.get('SAMPLE',''),
                            acquired_at=time.time())
                        _save(path,value);leases.append(Lease(device,path,value['reserved_bytes']))
                except BaseException:
                    for lease in leases:lease.path.unlink(missing_ok=True)
                    raise
                return leases[0] if count==1 else LeaseGroup(leases)
        if time.monotonic()>=deadline:raise TimeoutError('GPU/VRAM 예약 대기 시간 초과: '+purpose)
        time.sleep(1)

@contextmanager
def reserved(required_bytes, *, purpose, exclusive=False):
    lease=acquire(required_bytes,purpose=purpose,exclusive=exclusive)
    try:yield lease
    finally:lease.close()
