"""실제 두 HAP STAR index의 직렬/병렬 동등성과 시간을 기록함."""
import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path

CORE = ('Genome','SA','SAindex','chrName.txt','chrLength.txt','chrStart.txt',
        'sjdbInfo.txt','sjdbList.out.tab','sjdbList.fromGTF.out.tab',
        'exonGeTrInfo.tab','exonInfo.tab','geneInfo.tab','transcriptInfo.tab')


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as source:
        while block:=source.read(8<<20):h.update(block)
    return h.hexdigest()


def snapshot(work):
    result={}
    for h in [1,2]:
        path=work/f'hap{h}.index.json';record=json.loads(path.read_text());index=Path(record['index_dir'])
        hashes={name:{'bytes':(index/name).stat().st_size,'sha256':digest(index/name)}
                for name in CORE if (index/name).is_file()}
        if not all(name in hashes for name in CORE[:6]):raise ValueError('STAR index 핵심 파일 누락')
        result[f'HAP{h}']={'record':record,'hashes':hashes}
    return result


def save(path,record):
    temporary=path.with_name('.'+path.name+f'.{os.getpid()}.part')
    temporary.write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n');os.replace(temporary,path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['prepare','snapshot','compare']);parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True);parser.add_argument('--seconds',type=float,required=True)
    args=parser.parse_args()
    if args.mode=='prepare':
        # 두 경로 모두 target FASTA 구성 비용까지 동일하게 다시 측정함.
        for h in [1,2]:
            path=args.work_dir/f'hap{h}.target.json'
            if path.is_file():
                record=json.loads(path.read_text());record['cache_fingerprint']='benchmark_force_target'
                save(path,record)
        return
    current=snapshot(args.work_dir)
    if args.mode=='snapshot':
        report={'schema':'mfm_star_index_benchmark_v1','status':'serial_complete','serial_seconds':args.seconds,'serial':current,
                'scope':'동일한 개인 HAP·GTF·추가 contig로 생성한 핵심 STAR index 파일의 SHA256 동등성 비교임. 측정용 두 번째 생성 시간과 해시 시간을 일반 전처리 시간과 구분함.'}
    else:
        report=json.loads(args.report.read_text());same=all(report['serial'][h]['hashes']==current[h]['hashes'] for h in current)
        intervals=[current[h]['record'].get('index_execution',{}) for h in ['HAP1','HAP2']]
        serial_runs=[report['serial'][h]['record'].get('index_execution',{}) for h in ['HAP1','HAP2']]
        try:
            overlap=max(0.0,min(float(x['finished_at']) for x in intervals)-max(float(x['started_at']) for x in intervals))
            valid=all(float(x['finished_at'])>float(x['started_at']) and int(x['threads'])>0 for x in [*intervals,*serial_runs])
            comparable=valid and sum(int(x['threads']) for x in intervals)<=min(int(x['threads']) for x in serial_runs)
        except (KeyError,ValueError,TypeError):overlap=0.0;comparable=False
        eligible=same and overlap>0 and comparable
        report.update(status='pass' if same else 'serial_fallback',parallel_seconds=args.seconds,
                      parallel=current,core_sha256_equal=same,
                      actual_overlap_seconds=overlap,cpu_budget_comparable=comparable,
                      recommended_hap_jobs=2 if eligible and args.seconds<report['serial_seconds']*.95 else 1)
        # 검증이 실패하거나 더 느린 경우 검증된 직렬 캐시를 사용함.
        if report['recommended_hap_jobs']==1:
            for h in [1,2]:save(args.work_dir/f'hap{h}.index.json',report['serial'][f'HAP{h}']['record'])
        if args.seconds:report['speedup']=report['serial_seconds']/args.seconds
    args.report.parent.mkdir(parents=True,exist_ok=True)
    label='serial' if args.mode=='snapshot' else 'parallel'
    logs=args.report.parent/'rna_index_benchmark_logs'/label
    for source in (args.work_dir.parent/'logs').glob('*.index.HAP*.log'):
        logs.mkdir(parents=True,exist_ok=True);shutil.copy2(source,logs/source.name)
    save(args.report,report)
    print('[STAR-INDEX-BENCH] '+json.dumps({k:v for k,v in report.items() if k not in {'serial','parallel'}},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
