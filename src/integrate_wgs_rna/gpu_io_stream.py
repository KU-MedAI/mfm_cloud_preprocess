"""최대 두 배치만 미리 읽어 IO와 GPU 계산을 겹침. 실패 시 생산자를 중단함."""
from queue import Queue,Full
from threading import Thread,Event


def prefetch(iterator,buffers=2):
    queue=Queue(maxsize=buffers);stop=Event()
    def send(value):
        while not stop.is_set():
            try:queue.put(value,timeout=.1);return
            except Full:pass
    def produce():
        try:
            for value in iterator:
                if stop.is_set():break
                send((True,value))
        except BaseException as error:send((False,error))
        finally:send((True,None))
    thread=Thread(target=produce,daemon=True);thread.start()
    try:
        while True:
            success,value=queue.get()
            if not success:raise value
            if value is None:break
            yield value
    finally:
        stop.set();thread.join(timeout=2)
