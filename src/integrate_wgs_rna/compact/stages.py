"""compact 생성·최종 주석을 독립된 완료 단위로 관리함."""
from pathlib import Path
from ..pipeline.checkpoints import run_step,file_identity,tree_identity,digest


def compact_step(sample_root, compact_root, records, config, action):
    root=Path(compact_root)
    sources=sorted({str(r[k]) for r in records for k in
                    ('hap1_fasta','hap2_fasta','hap1_chain','hap2_chain','bam','hap1_alignment','hap2_alignment') if r.get(k)})
    for r in records:
        if r.get('window_index_dir'):sources.append(str(Path(r['window_index_dir'])/'windows.tsv'))
    # 정답/좌표뿐 아니라 CRAM 해독용 reference와 인덱스 변경도 재사용을 막음.
    sources += [str(Path(r['window_index_dir'])/'windows.json') for r in records
                if r.get('window_index_dir') and (Path(r['window_index_dir'])/'windows.json').is_file()]
    for source in list(sources):
        p=Path(source)
        for extra in [Path(str(p)+suffix) for suffix in ('.fai','.crai','.bai')]+[p.with_suffix('.crai'),p.with_suffix('.bai'),p.parent/'target.fa',p.parent/'target.fa.fai']:
            if extra.is_file():sources.append(str(extra))
    def inputs():return [file_identity(p) for p in sorted(set(sources))]
    def outputs():return {str(r['contig']):tree_identity(root/r['contig'],exclude=('windows.tsv',)) for r in records}
    module=Path(__file__).parents[1]
    code={str(p.relative_to(module)):digest(p.read_text()) for p in
          [module/'direct_compact_writer.py',module/'legacy_shared_hap_converter.py',module/'window_quality.py',
           module/'pipeline/checkpoints.py',Path(__file__),Path(__file__).with_name('packing.py')]}
    key={'export':config,'records':records,'code':code}
    return run_step(Path(sample_root)/'metadata/pipeline_steps/compact.json',key,action,outputs,inputs_identity=inputs)


def annotation_step(sample_root, compact_root, index_dir, records, contract):
    from ..window_policy import persist_window_annotations,validate_window_annotations
    if index_dir is None:return persist_window_annotations(sample_root,compact_root,index_dir,records,contract)
    import json
    from ..window_contract import ensure_window_n_annotations
    index=Path(index_dir)
    meta=json.loads((index/'windows.json').read_text()) if (index/'windows.json').exists() else {}
    if 'storage_policy' not in meta:
        return persist_window_annotations(sample_root,compact_root,index,records,contract)
    ensure_window_n_annotations(index/'windows.tsv')
    sources=[index/'windows.tsv',index/'windows.json']
    if (index/'windows_blacklist.tsv').is_file():sources.append(index/'windows_blacklist.tsv')
    sources += [Path(r[k]) for r in records for k in ('hap1_chain','hap2_chain') if r.get(k)]
    def inputs():return [file_identity(p) for p in sources]
    def outputs():
        check=validate_window_annotations(Path(sample_root))
        return {'validation':check,'policy':file_identity(Path(sample_root)/'metadata/window_policy.json')}
    def action():
        value=persist_window_annotations(sample_root,compact_root,index,records,contract)
        return {'status':'pass','annotations':value}
    module=Path(__file__).parents[1]
    key={'root':str(compact_root),'contract':contract,'policy_code':digest((module/'window_policy.py').read_text()),
         'n_contract_code':digest((module/'window_contract.py').read_text())}
    result=run_step(Path(sample_root)/'metadata/pipeline_steps/annotations.json',key,action,outputs,inputs_identity=inputs)
    return result['annotations']
