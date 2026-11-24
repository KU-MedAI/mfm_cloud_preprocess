"""SHA256 원값과 중복 개수를 보존하는 bounded 분할정렬 multiset."""
from collections import Counter
from pathlib import Path
import os,shutil,tempfile
import numpy as np


class PartitionedBag:
    """순차 spool 후 64-way 분할함. 큰 leaf는 다음 digest byte로 재분할함.

    기본 leaf 원자료64MiB, 불일치 집계 임시배열까지 최대8배를 예약함.
    digest를 XOR/축약하지 않으며 multiplicity까지 정렬한 원값으로 대조함.
    """
    def __init__(self,path,*,leaf_bytes=64<<20,batch_records=8192):
        if leaf_bytes<64 or batch_records<1:raise ValueError('multiset leaf/batch 한도 오류')
        self.root=Path(str(path)+'.partitioned');self.root.mkdir(exist_ok=False)
        self.leaf_bytes=int(leaf_bytes);self.batch_records=int(batch_records)
        self.total={'src':Counter(),'out':Counter()};self.pending={};self.handles={};self.files={};self.kinds={}
        self.stats={'backend':'partitioned_sha256','leaf_limit_bytes':self.leaf_bytes,'leaf_peak_raw_bytes':0,
                    'comparison_leaves':0,'partition_passes':0}

    def add(self,side,kind,digest):
        side='src' if side=='src' else 'out';digest=bytes(digest)
        if len(digest)!=32:raise ValueError('SHA256 digest는32바이트여야 함')
        if kind not in self.kinds:self.kinds[kind]=len(self.kinds)
        key=(side,kind)
        if key not in self.files:self.files[key]=self.root/f'{self.kinds[kind]}.{side}.bin'
        data=self.pending.setdefault(key,bytearray());data.extend(digest);self.total[side][kind]+=1
        if len(data)>=32*self.batch_records:self._flush(key)

    def add_many(self,side,kind,digests):
        values=np.ascontiguousarray(digests)
        if values.ndim!=1 or values.dtype!=np.dtype('V32'):raise ValueError('SHA256 batch 오류')
        side='src' if side=='src' else 'out';key=(side,kind)
        if kind not in self.kinds:self.kinds[kind]=len(self.kinds)
        if key not in self.files:self.files[key]=self.root/f'{self.kinds[kind]}.{side}.bin'
        self.pending.setdefault(key,bytearray())
        self._flush(key)
        if key not in self.handles:self.handles[key]=self.files[key].open('ab',buffering=1<<20)
        self.handles[key].write(values.tobytes());self.total[side][kind]+=len(values)

    def _flush(self,key):
        data=self.pending[key]
        if data:
            if key not in self.handles:self.handles[key]=self.files[key].open('ab',buffering=1<<20)
            self.handles[key].write(data);data.clear()

    def _seal(self):
        for key in self.pending:self._flush(key)
        for handle in self.handles.values():handle.close()
        self.handles.clear()

    @staticmethod
    def _size(path):return path.stat().st_size if path is not None else 0

    def _split(self,path,depth,directory,label):
        if path is None:return {}
        outputs={};handles={}
        try:
            with path.open('rb') as source:
                while raw:=source.read(32768*32):
                    if len(raw)%32:raise ValueError('잘린 SHA256 spool')
                    rows=np.frombuffer(raw,dtype='V32');codes=np.frombuffer(raw,dtype=np.uint8).reshape(-1,32)[:,depth]
                    buckets=codes>>2 if depth==0 else codes
                    order=np.argsort(buckets,kind='stable');counts=np.bincount(buckets,minlength=256);ordered=rows[order]
                    offset=0
                    for bucket in np.flatnonzero(counts):
                        n=int(counts[bucket]);bucket=int(bucket)
                        if bucket not in handles:
                            outputs[bucket]=directory/f'{label}.{bucket}.bin';handles[bucket]=outputs[bucket].open('wb',buffering=1<<16)
                        handles[bucket].write(ordered[offset:offset+n].tobytes());offset+=n
        finally:
            for handle in handles.values():handle.close()
        return outputs

    @staticmethod
    def _new_item():return {'identity_count':0,'source_only':0,'zarr_only':0,'examples':[]}

    @staticmethod
    def _mismatch(item,digest,src,out,max_examples):
        if src==out:return
        item['identity_count']+=1;item['source_only']+=max(0,src-out);item['zarr_only']+=max(0,out-src)
        item['examples'].append({'sha256':bytes(digest).hex(),'source_count':src,'zarr_count':out})
        # recursive shard의 순서와 무관하게 전역 digest 순으로 표본을 유지함.
        item['examples'].sort(key=lambda x:x['sha256'])
        del item['examples'][max(0,max_examples):]

    def _leaf(self,src,out,item,max_examples):
        size=self._size(src)+self._size(out)
        self.stats['leaf_peak_raw_bytes']=max(self.stats['leaf_peak_raw_bytes'],size)
        self.stats['comparison_leaves']+=1
        a=np.fromfile(src,dtype='V32') if src else np.empty(0,dtype='V32')
        b=np.fromfile(out,dtype='V32') if out else np.empty(0,dtype='V32')
        a.sort();b.sort()
        if np.array_equal(a,b):return
        ua,ca=np.unique(a,return_counts=True);ub,cb=np.unique(b,return_counts=True)
        keys=np.union1d(ua,ub);ac=np.zeros(len(keys),dtype=np.int64);bc=ac.copy()
        ac[np.searchsorted(keys,ua)]=ca;bc[np.searchsorted(keys,ub)]=cb
        for i in np.flatnonzero(ac!=bc):self._mismatch(item,keys[i],int(ac[i]),int(bc[i]),max_examples)

    def _terminal(self,src,out,item,max_examples):
        # first byte상위6bit와나머지31byte가같으면최대4개정체성만남음.
        counts=[]
        for path in [src,out]:
            counter=Counter()
            if path:
                with path.open('rb') as source:
                    while raw:=source.read(32768*32):
                        for value in np.frombuffer(raw,dtype='V32'):counter[bytes(value)]+=1
                        if len(counter)>4:raise ValueError('terminal digest partition 오류')
            counts.append(counter)
        for digest in sorted(counts[0].keys()|counts[1].keys()):
            self._mismatch(item,digest,counts[0][digest],counts[1][digest],max_examples)

    def _compare(self,src,out,item,max_examples,depth,directory,remove=False):
        size=self._size(src)+self._size(out)
        if self._size(src)%32 or self._size(out)%32:raise ValueError('잘린 SHA256 spool')
        if size<=self.leaf_bytes:self._leaf(src,out,item,max_examples)
        elif depth>=32:self._terminal(src,out,item,max_examples)
        else:
            self.stats['partition_passes']+=1
            node=directory/f'node{self.stats["partition_passes"]}';node.mkdir()
            left=self._split(src,depth,node,'src');right=self._split(out,depth,node,'out')
            # 중첩 분할 중 원래 spool과두단계 partition만 동시에 보존함.
            if remove:
                for p in [src,out]:
                    if p is not None:p.unlink()
                remove=False
            for bucket in sorted(left.keys()|right.keys()):
                self._compare(left.get(bucket),right.get(bucket),item,max_examples,depth+1,directory,True)
            node.rmdir()
        if remove:
            for p in [src,out]:
                if p is not None:p.unlink()

    def result(self,max_examples):
        self._seal();mismatch={};kinds=sorted(self.total['src'].keys()|self.total['out'].keys())
        with tempfile.TemporaryDirectory(prefix='compare_',dir=self.root) as tmp:
            for kind in kinds:
                item=self._new_item();self._compare(self.files.get(('src',kind)),self.files.get(('out',kind)),item,int(max_examples),0,Path(tmp))
                if item['identity_count']:mismatch[kind]=item
        return {'pass':not mismatch,'source_total':{k:int(self.total['src'][k]) for k in kinds},
                'zarr_total':{k:int(self.total['out'][k]) for k in kinds},'mismatches':mismatch}

    def close(self):
        self._seal()
        if self.root.exists():shutil.rmtree(self.root)
