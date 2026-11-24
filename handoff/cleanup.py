"""검증된 전달본에서 모델 입력·QC·메타데이터·로그만 보존한다."""
from __future__ import annotations
import argparse, json, os, shutil, tempfile, time
from pathlib import Path

from scripts.cleanup_guard import validate_predelete
from tools.compact_manifest import completed_contigs

KEEP_DIRS={'compact','metadata','qc','logs','rna_personal'}
KEEP_RNA={'coverage.zarr','coverage_windows.zarr'}
RAW_SUFFIXES=('.bam','.cram','.sam','.fa','.fasta','.fna','.fastq','.fq','.bam.bai','.cram.crai','.fa.fai','.fasta.fai','.fa.gz','.fasta.gz','.fastq.gz','.fq.gz','.fa.gz.fai','.fa.gz.gzi')

def model_files(root):
    root=Path(root)
    for base in ('compact','rna_personal'):
        for p in (root/base).rglob('*') if (root/base).is_dir() else []:
            if p.is_file() and not p.is_symlink():
                if base=='compact' or any(x in KEEP_RNA for x in p.parts) or p.suffix=='.tsv':yield p

def model_bytes(root):return sum(p.stat().st_size for p in model_files(root))

def _write_json_atomic(path, payload):
    """진행·실패·완료 보고서를 같은 디렉터리에서 원자적으로 기록한다."""
    path=Path(path)
    with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False,encoding='utf-8') as handle:
        json.dump(payload,handle,ensure_ascii=False,indent=2)
        handle.write('\n')
        temporary=Path(handle.name)
    try:
        os.replace(temporary,path)
    finally:
        if temporary.exists():temporary.unlink()

def cleanup_complete(root, sample, mode, *, external=True):
    """파일 존재가 아닌 완료 보고서와 정리 후 PASS 계약을 확인한다."""
    meta=Path(root)/'metadata'
    if (meta/'cleanup_in_progress.json').exists():return False
    if external and (meta/'external_cleanup_in_progress.json').exists():return False
    try:
        report=json.loads((meta/'cleanup_completed.json').read_text(encoding='utf-8'))
        if not isinstance(report,dict):return False
        post=report.get('post_cleanup_validation')
        return (report.get('schema')=='mfm_delivery_cleanup_v1'
                and report.get('sample')==sample and report.get('mode')==mode
                and report.get('dry_run') is False
                and report.get('status','complete')=='complete'
                and isinstance(post,dict)
                and post.get('schema')=='mfm_delivery_final_validation_v1'
                and post.get('sample')==sample and post.get('mode')==mode
                and post.get('status')=='pass' and post.get('errors')==[])
    except (OSError,ValueError):
        return False

def cleanup(root, sample, mode, work=None, output_root=None, dry_run=False):
    root=Path(root).absolute()
    if root.is_symlink() or root.name!=sample or root==Path('/'):raise ValueError('샘플 루트 경로 불일치')
    # 예전 PASS JSON만 신뢰하지 않고 삭제 직전 현재 RNA 근거의 schema/hash/freshness를 직접 본다.
    predelete_guard=validate_predelete(root,sample,mode)
    layouts=list((root/'compact').glob('*/manifest.json'))
    if len(layouts)!=1:raise ValueError('compact manifest는 하나여야 합니다.')
    m=json.loads(layouts[0].read_text())
    if m.get('status')!='complete' or m.get('validation_level')!='full':raise ValueError('full compact 검증 필요')
    traces=list((root/'qc/stage12_trace').glob('*/trace_manifest.json'))
    if len(traces)!=1:raise ValueError('독립 source trace 필요')
    t=json.loads(traces[0].read_text())
    expected=completed_contigs(m,layouts[0].parent,sample)
    if not expected or t.get('status')!='pass' or {k for k,v in t.get('contigs',{}).items() if v=='pass'}!=expected:
        raise ValueError('source trace contig 불일치')
    if Path(t.get('compact_root','')).resolve()!=layouts[0].parent.resolve():raise ValueError('source trace 경로 불일치')
    if any(p.is_file() and p.stat().st_mtime_ns>traces[0].stat().st_mtime_ns for p in layouts[0].parent.rglob('*')):
        raise ValueError('source trace 이후 compact가 변경되었습니다.')
    # symlink를 따라 외부 파일을 삭제하거나 보존 검사에서 빠뜨리지 않는다.
    if any((root/base).is_symlink() for base in KEEP_DIRS) or any(p.is_symlink() for base in KEEP_DIRS for p in (root/base).rglob('*') if (root/base).exists()):
        raise ValueError('최종 보존 파일에 외부 symlink가 있습니다.')
    meta=root/'metadata';meta.mkdir(exist_ok=True)
    remove=[]; snapshots=[]
    for p in sorted(root.iterdir()):
        if p.name not in KEEP_DIRS:remove.append(p)
    for p in (root/'rna_personal').rglob('*') if (root/'rna_personal').exists() else []:
        if p.is_file() and not any(x.endswith('.zarr') for x in p.parts):
            if p.suffix not in {'.json','.tsv','.log','.txt'}:remove.append(p)
    # BAM/FASTA는 QC/logs 하위에 잘못 들어간 경우에도 모델 전용 패키지에서 제거한다.
    for base in KEEP_DIRS:
        for p in (root/base).rglob('*') if (root/base).exists() else []:
            if p.is_file() and p.name.lower().endswith(RAW_SUFFIXES):remove.append(p)
    fasta=root/'individual_fasta'
    for name in ['permanent_assets.json','personal_fasta_qc.json']:
        if (fasta/name).is_file():snapshots.append((fasta/name,meta/name))
    if (root/'windows/windows.tsv').is_file():snapshots.append((root/'windows/windows.tsv',meta/'wgs_windows.tsv'))
    if work:
        given=Path(work).absolute();w=given.resolve();o=Path(output_root).resolve() if output_root else None
        if not o or given.is_symlink() or o not in w.parents or '_work' not in w.relative_to(o).parts:
            raise ValueError('work 경로가 전달본 output/_work 아래가 아닙니다.')
        if w.exists():remove.append(w)
        if (w/'logs').is_dir() and not dry_run:shutil.copytree(w/'logs',root/'logs/work',dirs_exist_ok=True)
    unique=sorted(set(remove),key=lambda p:len(p.parts))
    targets=[]
    for p in unique:
        if not any(parent in targets for parent in p.parents):targets.append(p)
    payload={'schema':'mfm_delivery_cleanup_v1','sample':sample,'mode':mode,'dry_run':dry_run,
      'time':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'policy':'모델 입력/QC/metadata/logs 보존; sample BAM/CRAM/FASTA/chain 제거',
      'remove':[str(p) for p in targets],'model_bytes':model_bytes(root),
      'predelete_guard':predelete_guard,
      'resume_note':'정리 후 RNA만 재실행할 수 없다. 원본 다운로드와 WGS부터 별도 실행해야 한다.'}
    if (meta/'sample.json').is_file() and json.loads((meta/'sample.json').read_text()).get('rna_processing_status')=='deferred':
        payload['completion_scope']='dna_pileup_feature_only_rna_deferred' if json.loads((meta/'sample.json').read_text()).get('wgs_payload')=='pileup_feature' else 'dna_sequence_only_rna_deferred'
        payload['resume_note']='DNA 서열 Zarr 완료, RNA/ASE/CNV 미처리임. 추후 동일 reference/VCF/window 좌표로 read 근거를 재구성하고 추가 검증해야 함.'
    if dry_run:
        payload['status']='planned'
        _write_json_atomic(meta/'cleanup_plan.json',payload)
        return payload
    progress=meta/'cleanup_in_progress.json'
    payload['status']='running'
    _write_json_atomic(progress,payload)
    try:
        for a,b in snapshots:shutil.copy2(a,b)
        for p in targets:
            if p.is_symlink() or p.is_file():p.unlink()
            elif p.exists():shutil.rmtree(p)
        from .final_check import validate
        payload['post_cleanup_validation']=validate(root,sample,mode)
        payload['status']='complete'
        _write_json_atomic(meta/'cleanup_completed.json',payload)
        progress.unlink()

    except Exception as error:
        payload['status']='failed'
        payload['error']=str(error)
        _write_json_atomic(progress,payload)
        raise
    if (os.environ.get('MFM_MODEL_ONLY_CLEANUP') == '1'
            and os.environ.get('INPUT_ROOT') and os.environ.get('OUTPUT_ROOT')):
        from .completion_cleanup import cleanup_external
        cleanup_external(root, sample, mode, os.environ['INPUT_ROOT'],
                         os.environ['OUTPUT_ROOT'],
                         allow_unowned=os.environ.get('MFM_EXTERNAL_CLEANUP_UNOWNED') == '1')
    return payload

if __name__=='__main__':
    p=argparse.ArgumentParser(description='최종 모델 입력 패키지 정리')
    for name in ['root','sample','mode','work','output-root']:p.add_argument('--'+name,required=name in ['root','sample','mode'])
    p.add_argument('--dry-run',action='store_true');a=p.parse_args()
    print(json.dumps(cleanup(a.root,a.sample,a.mode,a.work,a.output_root,a.dry_run),ensure_ascii=False))
