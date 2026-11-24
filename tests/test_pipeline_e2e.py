"""실제 BAM/CRAM·FASTA→compact→주석→독립trace→package→정리 E2E."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,tempfile,unittest
from unittest.mock import patch
import pysam
import zarr
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src'),str(ROOT/'vendor/MFM-DNARNAFM')]
from integrate_wgs_rna.full_preprocess import FullPreprocessConfig,run_full_preprocess
from integrate_wgs_rna.personal_fasta import reference_fingerprint
from integrate_wgs_rna.window_index import WindowIndexConfig,run_window_index_generation
from integrate_wgs_rna.window_policy import validate_window_annotations
from integrate_wgs_rna.stage70_trace_validator import validate_stage70_trace_root
from integrate_wgs_rna.model_input import load_window
from integrate_wgs_rna.window_reader import WindowReader


class PipelineE2ETest(unittest.TestCase):
    def test_real_engine_resume_trace_and_cleanup(self):
        with tempfile.TemporaryDirectory() as raw:
            base=Path(raw);sample=base/'TEST';sample.mkdir()
            sequence='NN'+'ACGT'*5+'N'+'ACGT'*5+'N'*32+'ACGT'*7+'NN'
            ref=base/'reference.fa';ref.write_text('>chr1\n'+sequence+'\n');pysam.faidx(str(ref))
            alignment=base/'source.bam'
            header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':len(sequence),'M5':hashlib.md5(sequence.encode()).hexdigest()}]}
            with pysam.AlignmentFile(str(alignment),'wb',header=header) as out:
                for i,pos in enumerate([3,7,22,78]):
                    read=pysam.AlignedSegment();read.query_name=f'r{i}';read.query_sequence='ACGT';read.query_qualities=pysam.qualitystring_to_array('IIII')
                    read.flag=0;read.reference_id=0;read.reference_start=pos;read.mapping_quality=60;read.cigar=[(0,4)];read.set_tag('HP',1+i%2);read.set_tag('PS',1);out.write(read)
            pysam.index(str(alignment))
            results={}
            for hap in ['HAP1','HAP2']:
                d=sample/'individual_fasta'/hap;d.mkdir(parents=True)
                fa=d/'chr1.fa';fa.write_text(ref.read_text());pysam.faidx(str(fa))
                (d/'chr1.chain').write_text(f'chain 0 chr1 {len(sequence)} + 0 {len(sequence)} chr1 {len(sequence)} + 0 {len(sequence)} 1\n{len(sequence)}\n')
                realigned=sample/'realigned'/hap;realigned.mkdir(parents=True)
                target=realigned/'target.fa';target.write_text(ref.read_text());pysam.faidx(str(target))
                cram=realigned/'source.cram'
                with pysam.AlignmentFile(str(alignment),'rb') as source, pysam.AlignmentFile(str(cram),'wc',header=source.header,reference_filename=str(target)) as out:
                    for read in source:out.write(read)
                pysam.index(str(cram));results[hap]=[{'contig':'chr1','chain_qc':{'status':'pass'}}]
            qc={'schema':'mfm_personal_fasta_sample_qc_v1','sample':'TEST','status':'pass','sequence_variant_policy':'snp_only',
                'reference_fasta':str(ref),'reference_fingerprint':reference_fingerprint(ref),'hap_results':results}
            (sample/'individual_fasta/personal_fasta_qc.json').write_text(json.dumps(qc))
            run_window_index_generation(WindowIndexConfig(cohort='TEST',output_root=base,samples=['TEST'],chroms=['chr1'],
                sequence_variant_policy='snp_only',window_size=16,offset_unit=1,blacklist_bed=None))
            config=FullPreprocessConfig(sample_id='TEST',dataset_role='pre_training',read_type='short_read',project_id='TEST',
                hap1_fasta_dir=sample/'individual_fasta/HAP1',hap1_chain_dir=sample/'individual_fasta/HAP1',
                hap2_fasta_dir=sample/'individual_fasta/HAP2',hap2_chain_dir=sample/'individual_fasta/HAP2',
                alignment=alignment,hap1_alignment=sample/'realigned/HAP1/source.cram',hap2_alignment=sample/'realigned/HAP2/source.cram',
                reference_fasta=ref,sequence_reference_fasta=ref,personal_fasta_qc=sample/'individual_fasta/personal_fasta_qc.json',
                sequence_variant_policy='snp_only',contigs=['chr1'],output_sample_root=sample,window_index_dir=sample/'windows',
                window_size=16,offset_unit=1,jobs=1,windows_per_chunk=2,layout_name='layout',storage_policy='keep_raw_intermediate')
            # 큰 계산 후 주석 단계만 실패해도 다음 호출이 실제 compact 엔진을 재실행하지 않음.
            with patch('integrate_wgs_rna.window_policy.persist_window_annotations',side_effect=RuntimeError('주석 실패 fixture')):
                with self.assertRaisesRegex(RuntimeError,'주석 실패'):run_full_preprocess(config)
            with patch('integrate_wgs_rna.full_preprocess.run_direct_compact_jobs',side_effect=AssertionError('이미 완료된 compact 재실행')):
                self.assertEqual(run_full_preprocess(config)['status'],'complete')
            self.assertEqual(validate_window_annotations(sample)['status'],'pass')
            layout=sample/'compact/layout';report=sample/'qc/stage12_trace/all_to_both'
            self.assertEqual(validate_stage70_trace_root(layout,['chr1'],report,workers=1)['summary']['status'],'pass')
            with patch('integrate_wgs_rna.stage70_trace_validator._validate_hap',side_effect=AssertionError('완료 HAP 재검증')):
                self.assertEqual(validate_stage70_trace_root(layout,['chr1'],report,workers=1)['summary']['status'],'pass')
            # upstream 정렬 단계는 준비된 fixture를 사용했음을 skip 이력으로 명시함.
            meta=sample/'metadata';logs=sample/'logs/stage_logs/WGS';logs.mkdir(parents=True,exist_ok=True)
            original=(ROOT/'09_personal_fasta.sh').exists();stop=13 if original else 11
            events=[{'schema':'mfm_execution_provenance_v1','sample_id':'TEST','branch':'WGS','stage':f'{i:02d}',
                     'event':'skip','status':'skipped','reason':'준비된 실제 BAM/FASTA fixture 사용'} for i in range(1,stop)]
            (meta/'execution_provenance.jsonl').write_text(''.join(json.dumps(x)+'\n' for x in events))
            (sample/'qc/wgs_preflight.json').write_text(json.dumps({'schema':'mfm_preflight_qc_v2','sample':'TEST','status':'pass'}))
            cmd=[sys.executable,str(ROOT/'scripts/validate_sample_package.py'),'--sample','TEST','--root',str(sample),'--mode','wgs',
                 '--layout','layout','--chroms','chr1','--hap-read-mode','all_to_both',
                 '--output',str(sample/'qc/final_package_validation.json'),'--completeness-tsv',str(meta/'completeness.tsv')]
            checked=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(checked.returncode,0,checked.stdout+checked.stderr)
            value=load_window(sample,'chr1','HAP1',1)
            self.assertFalse(value['dna_target_mask'][value['dna_codes']==4].any())
            from handoff.cleanup import cleanup
            cleanup(sample,'TEST','wgs',output_root=base)
            self.assertFalse((sample/'individual_fasta').exists())
            self.assertEqual(validate_window_annotations(sample)['status'],'pass')
            value=load_window(sample,'chr1','HAP1',1)
            self.assertEqual(len(value['wgs_depth']),16)
            reader=WindowReader(sample)
            cached=reader.load('chr1','HAP1',1)
            self.assertEqual(cached['dna_codes'].tolist(),value['dna_codes'].tolist())
            self.assertEqual(cached['wgs_depth'].tolist(),value['wgs_depth'].tolist())
            self.assertFalse(cached['dna_target_mask'][cached['dna_codes']==4].any())
            reader.close()


if __name__=='__main__':unittest.main()
