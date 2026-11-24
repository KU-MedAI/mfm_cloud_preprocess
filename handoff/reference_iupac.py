"""기존 모호염기→N 표현을 public VCF 입력 참조에도 명시적으로 맞춤."""
from pathlib import Path
import fcntl,hashlib,json,os

ALLOWED=b'ACGTNRYSWKMBDHV'
TABLE=bytes(ord('N') if i in b'RYSWKMBDHV' else i for i in range(256))

def signature(path):
 p=Path(path);s=p.stat()
 return dict(path=str(p.resolve()),size=s.st_size,mtime_ns=s.st_mtime_ns,ctime_ns=s.st_ctime_ns)

def verify_reference(target,rows):
 import pysam
 with pysam.FastaFile(str(target)) as fasta:
  if list(fasta.references)!=[row['contig'] for row in rows]:raise ValueError('정규화 참조 contig 불일치')
  for row in rows:
   if fasta.get_reference_length(row['contig'])!=row['length']:raise ValueError('정규화 참조 길이 불일치')
   digest=hashlib.sha256()
   for start in range(0,row['length'],4<<20):digest.update(fasta.fetch(row['contig'],start,min(row['length'],start+(4<<20))).encode())
   if digest.hexdigest()!=row['normalized_sequence_sha256']:raise ValueError('정규화 참조 내용 SHA 불일치')

def normalized_reference(source,target):
 import numpy as np,pysam
 source=Path(source);target=Path(target);target.parent.mkdir(parents=True,exist_ok=True)
 record_path=Path(str(target)+'.source.json')
 with Path(str(target)+'.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  identity=signature(source)
  if target.exists():
   record=json.loads(record_path.read_text())
   if (record.get('schema')!='mfm_reference_iupac_to_n_v1' or record.get('source')!=identity
       or record.get('target')!=signature(target) or record.get('index')!=signature(str(target)+'.fai')):
    raise ValueError('모호염기 정규화 참조 캐시의 원본/내용이 변경됨')
   verify_reference(target,record['contigs'])
   return str(target)
  temporary=Path(str(target)+'.part.'+str(os.getpid()));rows=[];block=60*65536
  try:
   with pysam.FastaFile(str(source)) as fasta,temporary.open('wb') as out:
    for chrom,length in zip(fasta.references,fasta.lengths):
     out.write(('>'+chrom+'\n').encode());before=hashlib.sha256();after=hashlib.sha256();counts=np.zeros(256,dtype=np.int64)
     for start in range(0,length,block):
      seq=fasta.fetch(chrom,start,min(length,start+block)).upper().encode('ascii')
      histogram=np.bincount(np.frombuffer(seq,dtype=np.uint8),minlength=256);counts+=histogram
      if any(histogram[i] for i in range(256) if i not in ALLOWED):raise ValueError('참조에 지원하지 않는 염기 문자가 있음')
      result=seq.translate(TABLE);before.update(seq);after.update(result)
      complete=len(result)//60
      wrapped=np.empty((complete,61),dtype=np.uint8);wrapped[:,:60]=np.frombuffer(result[:complete*60],dtype=np.uint8).reshape(complete,60);wrapped[:,60]=10
      out.write(wrapped.tobytes())
      if len(result)%60:out.write(result[complete*60:]+b'\n')
     rows.append(dict(contig=chrom,length=length,source_sequence_sha256=before.hexdigest(),normalized_sequence_sha256=after.hexdigest(),ambiguity_counts={chr(i):int(counts[i]) for i in b'RYSWKMBDHV' if counts[i]}))
   if signature(source)!=identity:raise ValueError('정규화 중 원본 참조가 변경됨')
   os.replace(temporary,target);pysam.faidx(str(target))
   verify_reference(target,rows)
   record=dict(schema='mfm_reference_iupac_to_n_v1',source=identity,target=signature(target),index=signature(str(target)+'.fai'),changed_bases=sum(sum(r['ambiguity_counts'].values()) for r in rows),canonical_acgtn_changes=0,contigs=rows)
   from .cli import atomic_json
   atomic_json(record_path,record)
   return str(target)
  finally:temporary.unlink(missing_ok=True)
