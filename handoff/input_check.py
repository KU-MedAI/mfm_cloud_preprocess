"""Stage 01 전 RNA 전체 해독을 확인하고 quickcheck 결과를 이후 검사에 재사용한다."""
from __future__ import annotations
import argparse, json, os, subprocess
from pathlib import Path

def fingerprint(path):
    p=Path(path).resolve();s=p.stat()
    if not p.is_file() or s.st_size==0: raise ValueError('입력 파일이 없거나 비었습니다: '+str(p))
    return {'path':str(p),'size':s.st_size,'mtime_ns':s.st_mtime_ns,'ctime_ns':s.st_ctime_ns,'inode':s.st_ino}

def check(path, output, sample, samtools):
    from .rna_integrity import validate_full
    full_report=os.environ.get('RNA_FULL_DECODE_REPORT') or str(output)+'.full.json'
    validate_full(path,full_report,sample,samtools,decode_reference=os.environ.get('RNA_DECODE_REFERENCE'))
    before=fingerprint(path)
    result=subprocess.run([samtools,'quickcheck','-v',str(path)],capture_output=True,text=True)
    if result.returncode: raise ValueError('RNA quickcheck 실패: '+result.stdout+result.stderr)
    if before!=fingerprint(path): raise ValueError('검사 중 입력 파일이 변경되었습니다.')
    payload={'schema':'mfm_stage01_quickcheck_v1','sample':sample,'stage':'01','status':'pass','input':before,'full_decode_report':str(full_report)}
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    return payload

def reuse(path, report, sample):
    d=json.loads(Path(report).read_text())
    if d.get('schema')!='mfm_stage01_quickcheck_v1' or d.get('status')!='pass' or d.get('sample')!=sample:
        raise ValueError('Stage 01 quickcheck 계약 불일치')
    if d.get('input')!=fingerprint(path): raise ValueError('Stage 01 이후 RNA 입력이 변경되었습니다.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description='Stage 01 RNA 입력 무결성 검사')
    for name in ['rna','output','sample','samtools']:p.add_argument('--'+name,required=True)
    a=p.parse_args();check(a.rna,a.output,a.sample,a.samtools)
