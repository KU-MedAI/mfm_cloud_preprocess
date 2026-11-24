"""독립 trace의 batch 4bit 해제와 기존 JSON 바이트 규약을 그대로 가속함."""
import hashlib,json
import numpy as np

_ASCII=np.full(16,ord('?'),dtype=np.uint8)
_ASCII[:5]=np.frombuffer(b'ACGTN',dtype=np.uint8)


def decode_reads(packed,offsets,lengths):
    packed=np.asarray(packed,dtype=np.uint8)
    codes=np.empty(packed.size*2,dtype=np.uint8)
    codes[::2]=packed>>4;codes[1::2]=packed&15
    text=_ASCII[codes].tobytes();base=int(offsets[0]);result=[]
    for i,length in enumerate(lengths):
        length=int(length);a=2*(int(offsets[i])-base);b=2*(int(offsets[i+1])-base)
        if a<0 or b<a or b>len(text) or length<0 or length>b-a:raise ValueError('invalid packed read')
        seq=text[a:a+length]
        if b'?' in seq:raise ValueError('invalid packed read')
        result.append(seq.decode('ascii'))
    return result


def read_digest(sequence,hp,ps,start,end,segment_type,cigar):
    # sequence는 위 decoder가 ACGTN만 반환하므로 JSON escape가 필요 없음.
    # 모든 숫자는 caller에서 기존 구현과 같이 Python int로 변환함.
    ops=','.join(map(str,cigar))
    encoded=f'["{sequence}",{hp},{ps},{start},{end},{segment_type},[{ops}]]'.encode('ascii')
    return hashlib.sha256(encoded).digest()


def segment_digest(read,index,q0,q1,h0,h1,strand):
    fields=(index,q0,q1,h0,h1,strand)
    if not all(type(x) is int for x in fields):
        encoded=json.dumps((bytes(read).hex(),*fields),separators=(',',':'),ensure_ascii=True).encode('utf-8')
    else:encoded=('["'+bytes(read).hex()+'",'+','.join(map(str,fields))+']').encode('ascii')
    return hashlib.sha256(encoded).digest()
