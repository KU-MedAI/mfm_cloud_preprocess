"""BAM/FASTA 정리 이후 모델 입력 패키지와 보존 파일의 구조를 다시 확인한다."""
from pathlib import Path
import json
from .cleanup import KEEP_DIRS,KEEP_RNA,package_inventory
from tools.compact_manifest import completed_contigs
from integrate_wgs_rna.validation_parallel import sequence_checks

def validate(root,sample,mode):
    import zarr
    root=Path(root);errors=[];arrays=0
    for p in root.iterdir():
        if p.name not in KEEP_DIRS:errors.append('불필요한 최상위 산출물: '+str(p))
    inventory=package_inventory(root)
    errors.extend('외부 symlink: '+str(p) for p in inventory['symlinks'])
    errors.extend('제거되지 않은 BAM/FASTA: '+str(p) for p in inventory['raw_files'])
    for p in inventory['zarr_dirs']:
        try:
            g=zarr.open_group(str(p),mode='r')
            # 데이터 내용은 정리 전 full/source trace가 검사했고 메타데이터는 모두 다시 읽음.
            counts=[]
            def visit(name,obj):
                if isinstance(obj,zarr.Array):counts.append(obj.size)
            g.visititems(visit)
            if not counts:errors.append('빈 Zarr: '+str(p))
            arrays+=len(counts)
        except Exception as e:errors.append(str(p)+': '+str(e))
    manifests=list((root/'compact').glob('*/manifest.json'))
    if len(manifests)!=1:errors.append('compact manifest 개수 오류')
    else:
        m=json.loads(manifests[0].read_text())
        try:expected=completed_contigs(m,manifests[0].parent,sample)
        except (OSError,ValueError,TypeError,AttributeError) as error:
            errors.append('compact 완료 기록 오류: '+str(error));expected=set()
        sequence_only=m.get('wgs_payload')=='sequence_only'
        pileup=m.get('wgs_payload')=='pileup_feature'
        ordered_contigs=sorted(expected)
        sequence_results=(sequence_checks([manifests[0].parent/c for c in ordered_contigs])
                          if sequence_only else [])
        if pileup:
            sequence_results=sequence_checks([manifests[0].parent/c for c in ordered_contigs])
        for index,contig in enumerate(ordered_contigs):
            d=manifests[0].parent/contig
            required=['hap_sequence.zarr','readset.HAP1.zarr','readset.HAP2.zarr','depth.zarr','windows.tsv','contig_metadata.json','filter_qc.json']
            if sequence_only:
                required=['hap_sequence.zarr','windows.tsv','contig_metadata.json','filter_qc.json']
                if sequence_results[index][1] is not None:
                    errors.append('서열전용 정리후 검증 실패: '+sequence_results[index][1])
            if pileup:
                required=['hap_sequence.zarr','pileup.zarr','windows.tsv','contig_metadata.json','filter_qc.json']
                if sequence_results[index][1] is not None:errors.append('pileup 정리후 검증 실패: '+sequence_results[index][1])
            for item in required:
                if not (d/item).exists():errors.append('필수 모델 입력 누락: '+str(d/item))
    if mode=='rna':
        labels=list((root/'rna_personal').glob('*'))
        if len(labels)!=1:errors.append('RNA evidence 디렉터리 개수 오류')
        else:
            for h in ['HAP1','HAP2']:
                for name in KEEP_RNA:
                    if not (labels[0]/h/name).is_dir():errors.append('RNA coverage 누락: '+h+'/'+name)
            if json.loads((root/'metadata/sample.json').read_text()).get('rna_payload','full')=='coverage_only':
                try:
                    from rna_signal import validate_signal
                    validate_signal(root,sample,labels[0].name,after_cleanup=True)
                except Exception as exc:errors.append('RNA signal 정리후 검증 실패: '+str(exc))
            else:
                for name in ['hap_variants.tsv','hap_variants_qc.json','allelic_counts.tsv','allelic_counts.both_haps.tsv','allelic_counts.split_half.tsv','allelic_counts.mmap_adjusted.tsv','ase_qc.json','ase_training_eligibility.tsv']:
                    if not (labels[0]/name).is_file():errors.append('RNA 최종 파일 누락: '+name)
    policy_path=root/'metadata/window_policy.json'
    candidate_meta=root/'metadata/wgs_windows.json'
    annotations=None
    try:
        new_policy=candidate_meta.exists() and json.loads(candidate_meta.read_text()).get('storage_policy')=='annotated_acgtn_v1'
        if policy_path.exists() or new_policy:
            from integrate_wgs_rna.window_policy import validate_window_annotations
            annotations=validate_window_annotations(root)
    except (OSError,ValueError,KeyError,TypeError) as error:
        errors.append('정리 후 window 주석·reference mapping 검증 실패: '+str(error))
    result={'schema':'mfm_delivery_final_validation_v1','sample':sample,'mode':mode,'status':'fail' if errors else 'pass',
      'errors':errors,'zarr_array_metadata_checked':arrays,'model_input_bytes':inventory['model_bytes'],
      'validation_scope':'정리 전 full+source trace PASS에 이어 정리 후 보존 구조/Zarr array 메타데이터 확인'}
    if annotations is not None:result['window_annotations_and_mapping']=annotations
    (root/'qc/delivery_final_validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    if errors:raise ValueError('정리 후 최종 검증 실패: '+'; '.join(errors[:5]))
    return result
