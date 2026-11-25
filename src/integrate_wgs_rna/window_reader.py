"""완료 패키지 전용 읽기 전용 reader. 메타데이터만 제한적으로 재사용함."""
from collections import OrderedDict
from pathlib import Path
import json
import os
import numpy as np
import zarr
from .model_input import dna_target_mask
from .window_chunk_cache import DecodedChunkCache


def _stamp(paths):
    result=[]
    for p in paths:
        s=p.stat();result.append((str(p),s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns))
    return tuple(result)


def _metadata_files(directory, arrays):
    paths=[directory/name for name in ('.zattrs','.zgroup') if (directory/name).is_file()]
    for name in arrays:
        paths.extend(p for p in (directory/name).iterdir() if p.is_file())
    return sorted(set(paths))


class WindowReader:
    """완료 후 불변인 sample에 사용함. fork/pickle 후에는 worker에서 새로 연결함.

    DNA/depth chunk는 설정한 바이트 상한 안에서만 재사용함. QC/manifest와 이미 읽은 좌표 파일의
    stat 지문 변경 시 실패함. 보관 파일 전체 내용 해시 검사와 구분함.
    """
    def __init__(self, sample_root, *, max_cached_contigs=2, max_metadata_bytes=64<<20, max_data_cache_bytes=64<<20):
        if max_cached_contigs<1 or max_metadata_bytes<1:raise ValueError('캐시 한도는 양수여야 함')
        self._config={'sample_root':str(Path(sample_root).resolve()),'max_cached_contigs':int(max_cached_contigs),
                      'max_metadata_bytes':int(max_metadata_bytes),'max_data_cache_bytes':int(max_data_cache_bytes)}
        self.root=Path(self._config['sample_root'])
        self._chunks=DecodedChunkCache(max_data_cache_bytes)
        self._pid=os.getpid();self._contigs=OrderedDict();self._rna={};self._bytes=0
        qc_path=self.root/'qc/final_package_validation.json'
        qc=json.loads(qc_path.read_text())
        if qc.get('status')!='pass' or qc.get('failed_checks') or qc.get('sample')!=self.root.name:
            raise ValueError('최종 검증을 통과한 모델 입력이 아님')
        manifests=list((self.root/'compact').glob('*/manifest.json'))
        if len(manifests)!=1:raise ValueError('모델 layout 개수 오류')
        manifest=json.loads(manifests[0].read_text())
        if manifest.get('status')!='complete' or manifest.get('validation_level')!='full':raise ValueError('미완료 compact')
        self.layout=manifests[0].parent
        self.sequence_only=manifest.get('wgs_payload')=='sequence_only'
        self.pileup=manifest.get('wgs_payload')=='pileup_feature'
        if self.pileup:
            from .pileup_features.schema import SCHEMA, LEGACY_SCHEMA, LEGACY_CHANNELS, resolve_channels
            self.pileup_schema=manifest.get('schema')
            self.pileup_feature_set=manifest.get('feature_set','legacy_reference12_v1')
            expected=(list(resolve_channels(manifest.get('channels'),self.pileup_feature_set))
                      if self.pileup_schema==SCHEMA else list(LEGACY_CHANNELS) if self.pileup_schema==LEGACY_SCHEMA else None)
            if expected is None or manifest.get('channels')!=expected:raise ValueError('pileup manifest schema/채널 불일치')
            self.pileup_channels=expected
        if manifest.get('wgs_payload','full') not in {'full','sequence_only','pileup_feature'}:raise ValueError('지원하지 않는 WGS payload')
        labels=[p for p in (self.root/'rna_personal').glob('*') if p.is_dir()]
        if len(labels)>1:raise ValueError('RNA label 개수 오류')
        self.rna_root=labels[0] if labels else None
        self._base_files=[qc_path,manifests[0]];self._base_stamp=_stamp(self._base_files)

    def __getstate__(self):return self._config
    def __setstate__(self,state):self.__init__(**state)

    def _check(self):
        if os.getpid()!=self._pid:self.__init__(**self._config)
        if _stamp(self._base_files)!=self._base_stamp:raise ValueError('완료 QC/manifest 변경: reader 재연결 필요')

    def close(self):
        self._contigs.clear();self._rna.clear();self._bytes=0;self._chunks.clear()

    def cache_info(self):
        return {'pid':self._pid,'cached_contigs':len(self._contigs),'cached_rna_haps':len(self._rna),
                'estimated_metadata_bytes':self._bytes,'decoded_cache_bytes':self._chunks.bytes,
                'chunk_hits':self._chunks.hits,'chunk_misses':self._chunks.misses,**{k:v for k,v in self._config.items() if k!='sample_root'}}

    def _reserve(self,size):
        while self._contigs and (len(self._contigs)>=self._config['max_cached_contigs'] or self._bytes+size>self._config['max_metadata_bytes']):
            _,entry=self._contigs.popitem(last=False);self._bytes-=entry['bytes']
        if self._bytes+size>self._config['max_metadata_bytes']:raise ValueError('window metadata 캐시 예산 초과')

    def _contig(self,contig):
        if contig in self._contigs:
            entry=self._contigs.pop(contig);self._contigs[contig]=entry
            if _stamp(entry['files'])!=entry['stamp']:raise ValueError('WGS 좌표 metadata 변경')
            return entry
        d=self.layout/contig
        sequence=zarr.open_group(str(d/'hap_sequence.zarr'),mode='r')
        depth=None if self.sequence_only else zarr.open_group(str(d/('pileup.zarr' if self.pileup else 'depth.zarr')),mode='r')
        if self.pileup:
            if depth.attrs.get('schema')!=self.pileup_schema or depth.attrs.get('channels')!=self.pileup_channels:raise ValueError('pileup schema/채널 불일치')
            from .pileup_features.schema import SCHEMA
            if self.pileup_schema==SCHEMA and depth.attrs.get('feature_set')!=self.pileup_feature_set:raise ValueError('pileup feature set 불일치')
        if sequence.attrs.get('sample')!=self.root.name:raise ValueError('DNA sample 불일치')
        estimated=sum(sequence[h]['lengths'].shape[0]*128+8192 for h in ('HAP1','HAP2'))
        self._reserve(estimated)
        entry={'hap':{},'bytes':estimated};files=_metadata_files(d/'hap_sequence.zarr',[])
        for hap in ('HAP1','HAP2'):
            seq=sequence[hap];dep=depth[hap] if depth is not None else None
            meta={k:seq[k][:] for k in ('offsets','lengths','actual_length','raw_start')}
            meta['depth_offsets']=dep['offsets'][:] if dep is not None else meta['offsets']
            n=len(meta['lengths'])
            if any(meta[k].shape!=(n,) for k in ('actual_length','raw_start')) or meta['offsets'].shape!=(n+1,) or meta['depth_offsets'].shape!=(n+1,):
                raise ValueError('WGS 좌표 metadata shape 불일치')
            if np.any(meta['actual_length']<=0) or np.any(meta['actual_length']>meta['lengths']) or np.any(meta['raw_start']<0):
                raise ValueError('WGS 실제 길이/시작점 오류')
            entry['hap'][hap]={'meta':meta,'bases':seq['bases_packed'],'depth':dep['depth'] if dep is not None else None}
            if self.pileup:entry['hap'][hap]['pileup']=dep
            files+=_metadata_files(d/'hap_sequence.zarr'/hap,('offsets','lengths','actual_length','raw_start'))
            files += [d/'hap_sequence.zarr'/hap/'bases_packed/.zarray']
            if dep is not None:
                store=d/('pileup.zarr' if self.pileup else 'depth.zarr')
                files+=_metadata_files(store/hap,('offsets',))
                files += [store/hap/'depth/.zarray']
                if self.pileup:
                    files+=_metadata_files(store,[])
                    from .pileup_features.schema import SCHEMA, ARRAYS
                    fields=ARRAYS if self.pileup_schema==SCHEMA else ('features','mapping_kind','baseq_depth','supplementary_depth')
                    files += [store/hap/name/'.zarray' for name in fields]
        entry['files']=files;entry['stamp']=_stamp(files)
        self._contigs[contig]=entry;self._bytes+=estimated
        return entry

    def _rna_hap(self,hap):
        if not self.rna_root:return None
        if hap in self._rna:
            entry=self._rna[hap]
            if _stamp(entry['files'])!=entry['stamp']:raise ValueError('RNA 좌표 metadata 변경')
            return entry
        d=self.rna_root/hap/'coverage_windows.zarr';group=zarr.open_group(str(d),mode='r')
        if group.attrs.get('sample')!=self.root.name or group.attrs.get('hap')!=hap:raise ValueError('RNA sample/HAP 불일치')
        n=group['contig'].shape[0];estimated=n*512+8192
        # RNA 상주 인덱스 때문에 기존 contig를 비울 수 있으므로 load는 필요한 참조를 로컬로 유지함.
        self._reserve(estimated)
        chrom=group['contig'][:];starts=group['personal_start'][:];ends=group['personal_end'][:];offsets=group['offsets'][:]
        if starts.shape!=(n,) or ends.shape!=(n,) or offsets.shape!=(n+1,):raise ValueError('RNA 좌표 shape 불일치')
        keys={}
        for i,(c,a,b) in enumerate(zip(chrom,starts,ends)):
            key=(str(c),int(a),int(b))
            if key in keys:raise ValueError('중복 RNA/WGS 대응 좌표')
            if int(b)<=int(a):raise ValueError('RNA 실제 길이 오류')
            keys[key]=i
        files=_metadata_files(d,('contig','personal_start','personal_end','offsets'))+[d/'depth/.zarray']
        entry={'keys':keys,'offsets':offsets,'depth':group['depth'],'files':files,'stamp':_stamp(files),'bytes':estimated}
        self._rna[hap]=entry;self._bytes+=estimated
        return entry

    def window_count(self,contig,hap):
        self._check();self._validate_key(contig,hap,0)
        return len(self._contig(contig)['hap'][hap]['meta']['lengths'])

    @staticmethod
    def _validate_key(contig,hap,index):
        if hap not in {'HAP1','HAP2'} or index<0:raise ValueError('HAP/window index 오류')
        if contig not in {f'chr{i}' for i in range(1,23)}:raise ValueError('상염색체 이름 오류')

    def load(self,contig,hap,index):
        self._check();self._validate_key(contig,hap,index)
        rna=self._rna_hap(hap)
        entry=self._contig(contig)['hap'][hap];meta=entry['meta']
        if index>=len(meta['lengths']):raise IndexError(index)
        begin,end=int(meta['offsets'][index]),int(meta['offsets'][index+1])
        padded=int(meta['lengths'][index]);actual=int(meta['actual_length'][index]);start=int(meta['raw_start'][index])
        if begin<0 or end-begin!=(padded+1)//2:raise ValueError('DNA packed offset 불일치')
        packed=self._chunks.read(entry['bases'],begin,end)
        if packed.size!=end-begin:raise ValueError('DNA packed 길이 불일치')
        codes=np.empty(packed.size*2,dtype=np.uint8);codes[0::2]=packed>>4;codes[1::2]=packed&15;codes=codes[:padded]
        mask=dna_target_mask(codes,actual)
        a,b=int(meta['depth_offsets'][index]),int(meta['depth_offsets'][index+1])
        values=None
        if entry['depth'] is not None:
            values=self._chunks.read(entry['depth'],a,b)
            if a<0 or b-a!=padded or values.shape!=codes.shape or values.dtype!=np.dtype('uint32'):raise ValueError('WGS 깊이 축/형식 불일치')
        result={'dna_codes':codes,'dna_target_mask':mask,'observation_mask':np.arange(padded)<actual,
                'wgs_depth':values,'actual_length':actual,'raw_start':start,'rna_depth':None}
        if self.pileup:
            g=entry['pileup']
            features=g['features'][a:b]
            kind=self._chunks.read(g['mapping_kind'],a,b)
            quality=self._chunks.read(g['baseq_depth'],a,b)
            if features.shape!=(padded,len(self.pileup_channels)) or features.dtype!=np.dtype('float32') or not np.isfinite(features).all():raise ValueError('pileup feature shape/dtype/값 오류')
            result.update(pileup_features=features,pileup_channels=list(self.pileup_channels),pileup_feature_set=self.pileup_feature_set,
                          pileup_reference_mapped_mask=kind==1,pileup_insertion_mask=kind==2,
                          pileup_baseq_valid_mask=quality>0,
                          supplementary_depth=self._chunks.read(g['supplementary_depth'],a,b),
                          max_deletion_length=self._chunks.read(g['max_deletion_length'],a,b),
                          wgs_depth_coordinate_system='reference_projected_to_hap')
            from .pileup_features.schema import SCHEMA
            if self.pileup_schema==SCHEMA:
                for name in ('hp1_depth','hp2_depth','hp1_supplementary_depth','hp2_supplementary_depth'):
                    result[name]=self._chunks.read(g[name],a,b)
                result.update(pileup_reference_base_code=self._chunks.read(g['reference_base_code'],a,b),
                              pileup_baf_valid_mask=self._chunks.read(g['baf_valid'],a,b)>0,
                              pileup_baf_site_mask=self._chunks.read(g['baf_site'],a,b)>0,
                              pileup_baf_heterozygous_mask=self._chunks.read(g['baf_heterozygous'],a,b)>0)
        if rna:
            try:i=rna['keys'][(contig,start,start+actual)]
            except KeyError as exc:raise ValueError('RNA/WGS 좌표 대응 오류') from exc
            a,b=int(rna['offsets'][i]),int(rna['offsets'][i+1]);observed=self._chunks.read(rna['depth'],a,b)
            if a<0 or b-a!=actual or observed.shape!=(actual,):raise ValueError('RNA 실제 길이 불일치')
            aligned=np.zeros(padded,dtype=observed.dtype);aligned[:actual]=observed;result['rna_depth']=aligned
        return result
