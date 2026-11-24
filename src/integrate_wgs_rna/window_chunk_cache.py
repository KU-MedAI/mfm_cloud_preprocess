"""1차원 Zarr의 해제된 chunk를 바이트 상한 안에서 재사용함."""
from collections import OrderedDict
from pathlib import Path
import numpy as np


def _stamp(path):
    try:s=path.stat()
    except FileNotFoundError:return None
    return s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns


class DecodedChunkCache:
    def __init__(self,limit):
        if limit<0:raise ValueError('chunk 캐시 상한은 음수가 될 수 없음')
        self.limit=int(limit);self.bytes=0;self.entries=OrderedDict();self.hits=0;self.misses=0

    def clear(self):self.entries.clear();self.bytes=0

    def read(self,array,start,end):
        if start<0 or end<start or end>array.shape[0]:raise ValueError('배열 조회 범위 오류')
        chunk=int(array.chunks[0]);chunk_bytes=chunk*array.dtype.itemsize
        if not self.limit or chunk_bytes>self.limit:return array[start:end]
        output=np.empty(end-start,dtype=array.dtype)
        for index in range(start//chunk,(end+chunk-1)//chunk):
            path=Path(array.store.path)/array.path/str(index)
            key=str(path);stamp=_stamp(path)
            if key in self.entries:
                data,old=self.entries.pop(key)
                if old!=stamp:
                    self.bytes-=data.nbytes
                    raise ValueError('캐시된 DNA/depth chunk 변경')
                self.hits+=1;self.entries[key]=(data,old)
            else:
                data=array.blocks[index].copy();self.misses+=1
                if stamp!=_stamp(path):raise ValueError('chunk를 읽는 동안 입력 변경')
                while self.entries and self.bytes+data.nbytes>self.limit:
                    _,(old,_) = self.entries.popitem(last=False);self.bytes-=old.nbytes
                self.entries[key]=(data,stamp);self.bytes+=data.nbytes
            a=max(start,index*chunk);b=min(end,index*chunk+len(data))
            output[a-start:b-start]=data[a-index*chunk:b-index*chunk]
        return output
