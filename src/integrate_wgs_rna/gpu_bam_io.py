"""samtools 정렬 comparator를 유지하고 최종 BGZF 압축을 GPU 자식에 연결함."""
from pathlib import Path
import json,os,subprocess,sys


def enabled():
    value=os.environ.get('MFM_IO_BACKEND','cpu')
    if value not in {'cpu','nvcomp'}:raise ValueError('MFM_IO_BACKEND는cpu/nvcomp여야함')
    return value=='nvcomp'


def sort_stream(view_command,output,*,samtools,sort_threads,sort_memory,log):
    output=Path(output);metadata=output.with_name(output.name+'.gpu_io.json')
    command=[samtools,'sort','-@',str(sort_threads),'-m',sort_memory,'-T',str(output)+'.sort',
             '--output-fmt','BAM,level=0','-']
    encoder=[sys.executable,'-B','-m','integrate_wgs_rna.gpu_bgzf','--output',str(output),
             '--metadata',str(metadata),'--verify-workers','1']
    processes=[]
    try:
        viewer=subprocess.Popen(view_command,stdout=subprocess.PIPE,stderr=log);processes.append(viewer)
        sorter=subprocess.Popen(command,stdin=viewer.stdout,stdout=subprocess.PIPE,stderr=log);processes.append(sorter)
        viewer.stdout.close()
        writer=subprocess.Popen(encoder,stdin=sorter.stdout,stdout=subprocess.DEVNULL,stderr=log);processes.append(writer)
        sorter.stdout.close()
        codes=[p.wait() for p in reversed(processes)]
        if any(codes):raise RuntimeError('GPU BAM IO 실패: writer/sort/view='+str(codes))
        result=json.loads(metadata.read_text())
        if not result.get('cpu_full_decode_and_byte_compare') or not result.get('cpu_crc_all_chunks_equal'):
            raise ValueError('GPU BAM IO의 CPU 전수 검증 근거 없음')
        return result
    except BaseException:
        for process in processes:
            if process.poll() is None:process.terminate()
        for process in processes:process.wait()
        raise
