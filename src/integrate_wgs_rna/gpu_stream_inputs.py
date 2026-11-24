"""검증된 FASTQ 캐시를 만들면서 bounded FIFO로 같은 바이트를 정렬기에 공급함."""
from pathlib import Path
from queue import Queue,Full,Empty
from threading import Thread,Event
import collections,errno,os,select,signal,tempfile,time
from .gpu_fastq import cached_inputs,MAX_LINE_BYTES

class _Replay:
    def __init__(self,lines,source):self.lines=collections.deque(lines);self.source=source
    def readline(self,size=-1):
        if not self.lines:return self.source.readline(size)
        line=self.lines.popleft()
        if size>=0 and len(line)>size:
            self.lines.appendleft(line[size:]);return line[:size]
        return line
    def read(self,size):
        if size<0:raise ValueError('Replay는 bounded read만 지원함')
        data=bytearray()
        while self.lines and len(data)<size:
            line=self.lines.popleft();n=min(len(line),size-len(data));data.extend(line[:n])
            if n<len(line):self.lines.appendleft(line[n:])
        if len(data)<size:data.extend(self.source.read(size-len(data)))
        return bytes(data)

class _Mirror:
    def __init__(self,path,stop):
        self.path=path;self.stop=stop;self.queue=Queue(2);self.error=None
        self.thread=Thread(target=self._run,daemon=True);self.thread.start()
    def _run(self):
        fd=None
        try:
            while fd is None and not self.stop.is_set():
                try:fd=os.open(self.path,os.O_WRONLY|os.O_NONBLOCK)
                except OSError as e:
                    if e.errno!=errno.ENXIO:raise
                    self.stop.wait(.05)
            if fd is None:return
            while not self.stop.is_set():
                try:data=self.queue.get(timeout=.1)
                except Empty:continue
                if data is None:return
                view=memoryview(data)
                while view and not self.stop.is_set():
                    try:view=view[os.write(fd,view):]
                    except BlockingIOError:select.select([], [fd], [], .1)
        except BaseException as e:self.error=e;self.stop.set()
        finally:
            if fd is not None:os.close(fd)
    def write(self,data):
        while not self.stop.is_set():
            try:self.queue.put(bytes(data) if data is not None else None,timeout=.1);return
            except Full:pass
        if self.error:raise self.error
        raise RuntimeError('FASTQ FIFO 공급 중단')
    def finish(self):
        self.write(None);self.thread.join()
        if self.error:raise self.error

def stream_cpu_budget(total,paired,producer_cpus):
    minimum=producer_cpus+(2 if paired else 1)+3
    requested=int(os.environ.get('MFM_GPU_STREAM_CPUS',str(max(minimum,min(16,total//4)))))
    if requested<minimum:raise ValueError('FIFO 생산자 최소 CPU보다 작은 설정임')
    return requested

class StreamInputs:
    def __init__(self,source,key,work,paired,budget,producer_cpus,expected_bytes=0):
        self.source=source;self.key=key;self.work=work;self.paired=paired;self.budget=budget
        self.producer_cpus=producer_cpus;self.expected_bytes=expected_bytes
        self.cpus=stream_cpu_budget(budget.cpus,paired,producer_cpus)
        if self.cpus>budget.cpus-2:raise ValueError('FIFO 생산자와 정렬 CPU 예산이 부족함')
        self.stop=Event();self.ready=Event();self.result=None;self.error=None;self.proc=None;self.mirrors=[];self.lengths=[]
        self.directory=tempfile.TemporaryDirectory(prefix='mfm_fifo_',dir=work)
        self.fifos=[Path(self.directory.name)/f'input_R{i+1}.fastq' for i in range(2 if paired else 1)]
        for path in self.fifos:os.mkfifo(path,0o600)
    def attach(self,proc):
        self.proc=proc;lines=[];size=0
        # profile 때문에 첫 read가 소실되지 않도록 같은 행을 다시 전달함.
        for _ in range(1000):
            first=proc.stdout.readline(MAX_LINE_BYTES+1)
            if not first:break
            record=[first]+[proc.stdout.readline(MAX_LINE_BYTES+1) for _ in range(3)]
            if any(len(v)>MAX_LINE_BYTES for v in record):raise MemoryError('FIFO profile 행 상한 초과')
            if not first.startswith(b'@') or not record[2].startswith(b'+') or not record[3] or len(record[1].rstrip(b'\r\n'))!=len(record[3].rstrip(b'\r\n')):raise ValueError('FIFO 최초 FASTQ 레코드 손상')
            lines.extend(record);size+=sum(map(len,record));self.lengths.append(len(record[1].rstrip(b'\r\n')))
            if size>=32*(1<<20):break
        if not self.lengths:raise ValueError('FIFO 입력 FASTQ가 비어 있음')
        self.mirrors=[_Mirror(path,self.stop) for path in self.fifos]
        self.ready.set();return _Replay(lines,proc.stdout)
    def finish_mirrors(self):
        for mirror in self.mirrors:mirror.finish()
    def _run(self):
        try:self.result=cached_inputs(self.source,self.key,self.work,self.paired,self.budget,self.producer_cpus,self.expected_bytes,streaming=self)
        except BaseException as e:self.error=e;self.stop.set()
        finally:self.ready.set()
    def __enter__(self):
        self.thread=Thread(target=self._run,daemon=True);self.thread.start()
        deadline=time.monotonic()+float(os.environ.get('MFM_GPU_WAIT_SECONDS','86400'))
        while not self.ready.wait(.1):
            if time.monotonic()>deadline:self.close();raise TimeoutError('FASTQ 최초 공급 준비 초과')
        if self.error:self.close();raise self.error
        self.reused=bool(self.result and self.result[2]);self.paths=self.result[0] if self.reused else self.fifos
        return self
    def check(self):
        if self.error:raise self.error
    def finish(self):
        self.thread.join()
        if self.error:raise self.error
        if self.result is None:raise RuntimeError('FIFO 완료 장부 없음')
        return self.result
    def close(self):
        self.stop.set()
        if self.proc is not None and self.proc.poll() is None:
            try:os.killpg(self.proc.pid,signal.SIGTERM)
            except ProcessLookupError:pass
        if hasattr(self,'thread'):self.thread.join(timeout=5)
        for mirror in self.mirrors:mirror.thread.join(timeout=2)
        if hasattr(self,'thread') and self.thread.is_alive():raise RuntimeError('FASTQ 생산자 종료 실패')
        self.directory.cleanup()
    def __exit__(self,kind,value,tb):
        try:
            if kind is None:self.finish()
        finally:self.close()
