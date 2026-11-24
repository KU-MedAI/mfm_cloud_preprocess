"""짧게 실행되는 GPU IO 검증 자식으로 GPU context/예약을 작업 종료 때 반환함."""
from pathlib import Path
import argparse,json,sys,time,os


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['cram-crc']);parser.add_argument('--path',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True);args=parser.parse_args()
    os.environ['MFM_IO_BACKEND']='nvcomp'
    started=time.monotonic()
    try:
        from .gpu_arrays import session,execution_metadata
        from .gpu_cram_integrity import scan
        with session('io') as cp:
            result=scan(args.path,cp)
            result['compute_execution']=execution_metadata()
        result.update(status='pass',elapsed_seconds=time.monotonic()-started)
    except Exception as error:
        from handoff.alignment_integrity import AlignmentCorruption
        result={'status':'failed','error':str(error),'corruption':isinstance(error,AlignmentCorruption),'error_type':type(error).__name__,
                'error_offset':getattr(error,'offset',None),'error_bytes':getattr(error,'length',None)}
    args.report.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':sys.exit(main())
