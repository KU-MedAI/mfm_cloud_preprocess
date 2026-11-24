"""실제 reference BAM과 consensus chain으로 수치·좌표·독립검증·정리 계약을 확인함."""
from pathlib import Path
from unittest.mock import patch
import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import pysam
import zarr

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src'),str(ROOT/'vendor/MFM-DNARNAFM')]
from integrate_wgs_rna.pileup_features.compute import calculate, ReferenceCache, projected
from integrate_wgs_rna.pileup_features.legacy_export import export, validate_contig, filter_config
from integrate_wgs_rna.pileup_features.schema import LEGACY_LAYOUT as LAYOUT
from integrate_wgs_rna.personal_fasta import reference_fingerprint
from integrate_wgs_rna.window_index import WindowIndexConfig,run_window_index_generation
from integrate_wgs_rna.window_reader import WindowReader
from tools.personal_assets import build,MANIFEST


def bam_file(path,sequence,records):
    header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':len(sequence),'M5':hashlib.md5(sequence.encode()).hexdigest()}]}
    with pysam.AlignmentFile(str(path),'wb',header=header) as out:
        for i,(seq,cigar,flag,mq,bq) in enumerate(records):
            r=pysam.AlignedSegment();r.query_name='r'+str(i);r.query_sequence=seq;r.flag=flag;r.reference_id=0;r.reference_start=0
            r.mapping_quality=mq;r.cigartuples=cigar
            if bq is not None:r.query_qualities=[bq]*len(seq)
            r.set_tag('HP',1+i%2);out.write(r)
    pysam.index(str(path))


class CountsTest(unittest.TestCase):
    def test_read_filter_config_matches_readset_defaults(self):
        from types import SimpleNamespace
        from integrate_wgs_rna.full_preprocess import stage70_export_config
        for read_type in ('short_read', 'long_read'):
            for min_mapq in (0, 30, 40):
                with self.subTest(read_type=read_type, min_mapq=min_mapq):
                    readset = stage70_export_config(SimpleNamespace(
                        read_type=read_type, stage70_filter_preset='auto',
                        window_size=1048576, offset_unit=None, seed=42, min_mapq=min_mapq))
                    self.assertEqual(filter_config(read_type, min_mapq).to_dict(), readset.to_dict())

    def test_clipping_padding_leading_indels_and_invalid_bases(self):
        with tempfile.TemporaryDirectory() as raw:
            path=Path(raw)/'edge.bam'
            records=[('TTACGT',[(5,2),(4,2),(0,2),(6,1),(0,2),(5,1)],16,60,20),
                     ('TTACGT',[(1,2),(0,4)],0,60,30),('ACGT',[(2,2),(0,4)],0,60,None)]
            bam_file(path,'ACGTACGT',records)
            with pysam.AlignmentFile(str(path),'rb') as bam:
                values=calculate(bam,'chr1',0,8,filter_config('short_read',30))
                oracle=calculate(bam,'chr1',0,8,filter_config('short_read',30),oracle=True)
                np.testing.assert_array_equal(values['depth'],[2,2,3,3,1,1,0,0])
                self.assertEqual(float(values['features'][:,9].sum()),0)
                for name in values:np.testing.assert_array_equal(values[name],oracle[name])
            bad=Path(raw)/'bad.bam';bam_file(bad,'ACGTACGT',[('ACRT',[(0,4)],0,60,30)])
            with pysam.AlignmentFile(str(bad),'rb') as bam:
                with self.assertRaisesRegex(ValueError,'염기'):calculate(bam,'chr1',0,8,filter_config('short_read',30))

    def test_mode_requires_reference_reads_and_separates_old_payloads(self):
        from handoff.sequence_mode import output_mode,selected,needs_wgs,runtime_environment,payload_matches,completion_mode
        from handoff.resume import sample_downloads
        row=dict(sample_id='S',read_type='short_read',data_group='wgs',alignment_path='wgs/S.bam',rna_path='')
        config=dict(data_root='/tmp/unused',pipeline_env={'MFM_SHORT_READ_OUTPUT_MODE':'pileup_feature'})
        self.assertEqual(output_mode(config,row),'pileup_feature');self.assertFalse(selected(config,row));self.assertTrue(needs_wgs(config,row))
        self.assertEqual(runtime_environment(config,row)['MFM_WGS_PAYLOAD'],'pileup_feature')
        resources=[dict(kind='wgs_alignment',path='wgs/S.bam')]
        self.assertEqual(sample_downloads(config,row,resources),resources)
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'qc').mkdir();q=root/'qc/final_package_validation.json'
            q.write_text(json.dumps({'wgs_payload':'sequence_only'}));self.assertFalse(payload_matches(config,row,root))
            q.write_text(json.dumps({'wgs_payload':'pileup_feature'}));self.assertFalse(payload_matches(config,row,root))
            from integrate_wgs_rna.pileup_features.schema import FEATURE_SET
            q.write_text(json.dumps({'wgs_payload':'pileup_feature','pileup_feature_set':FEATURE_SET}));self.assertTrue(payload_matches(config,row,root))
            pair=dict(row,data_group='wgs_rna_pair',rna_path='rna/S.bam')
            config['pipeline_env']['MFM_SEQUENCE_DEFER_RNA']='1'
            self.assertEqual(completion_mode(config,pair),'wgs');self.assertTrue(needs_wgs(config,pair))

    def test_indels_filters_quality_and_block_boundaries(self):
        with tempfile.TemporaryDirectory() as raw:
            path=Path(raw)/'source.bam'
            records=[('ACGT',[(0,4)],0,60,40),('ACTGT',[(0,2),(1,1),(0,2)],16,60,30),
                     ('AC',[(0,2),(2,2)],0,60,None),('ACNT',[(0,4)],0,60,0),
                     ('ACAC',[(0,2),(3,2),(0,2)],0,60,20)]
            records += [('ACGT',[(0,4)],flag,mq,40) for flag,mq in [(1024,60),(2048,60),(0,29),(256,60),(512,60)]]
            bam_file(path,'ACGTACGT',records)
            with pysam.AlignmentFile(str(path),'rb') as bam:
                config=filter_config('short_read',30);v=calculate(bam,'chr1',0,8,config)
                np.testing.assert_array_equal(v['depth'],[5,5,3,3,1,1,0,0])
                self.assertAlmostEqual(float(v['features'][1,9]),.2)
                self.assertAlmostEqual(float(v['features'][1,10]),.2)
                self.assertEqual(float(v['features'][1,11]),2)
                self.assertEqual(float(v['features'][0,2]),22.5)
                self.assertEqual(int(v['baseq_depth'][0]),4)
                self.assertAlmostEqual(float(v['features'][2,8]),1/3,places=6)
                for oracle in [False,True]:
                    other=calculate(bam,'chr1',0,8,config,oracle=oracle)
                    left=calculate(bam,'chr1',0,3,config,oracle=oracle);right=calculate(bam,'chr1',3,8,config,oracle=oracle)
                    for k in v:
                        np.testing.assert_array_equal(v[k],other[k])
                        np.testing.assert_array_equal(v[k],np.concatenate([left[k],right[k]],axis=0))
                long=calculate(bam,'chr1',0,8,filter_config('long_read',30))
                np.testing.assert_array_equal(long['depth'],[6,6,4,4,1,1,0,0])
                np.testing.assert_array_equal(long['supplementary_depth'],[1,1,1,1,0,0,0,0])


class PileupPackageTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name);self.root=self.base/'S1'
        self.sequence='NN'+'ACGT'*26+'N'+'ACGT'*12+'NN'
        self.ref=self.base/'reference.fa';self.ref.write_text('>chr1\n'+self.sequence+'\n');pysam.faidx(str(self.ref))
        n=len(self.sequence);self.hap2=self.sequence[:6]+'AT'+self.sequence[6:20]+self.sequence[22:]
        for hap,sequence in [('HAP1',self.sequence),('HAP2',self.hap2)]:
            fa=self.root/'individual_fasta'/hap/'chr1.fa';fa.parent.mkdir(parents=True);fa.write_text('>chr1\n'+sequence+'\n');pysam.faidx(str(fa))
            blocks=f'{n}\n' if hap=='HAP1' else f'6 0 2\n14 2 0\n{n-22}\n'
            fa.with_suffix('.chain').write_text(f'chain 0 chr1 {n} + 0 {n} chr1 {len(sequence)} + 0 {len(sequence)} 1\n'+blocks)
        qc=dict(schema='mfm_personal_fasta_sample_qc_v1',sample='S1',status='pass',sequence_variant_policy='snv_indel',
                reference_fasta=str(self.ref),reference_fingerprint=reference_fingerprint(self.ref),hap_results={})
        assets=self.root/'individual_fasta';(assets/'personal_fasta_qc.json').write_text(json.dumps(qc))
        m=build(assets,['chr1'],self.ref);self.assertEqual(m['status'],'pass');(assets/MANIFEST).write_text(json.dumps(m))
        bed=self.base/'blacklist.bed';bed.write_text('chr1\t2\t12\n')
        run_window_index_generation(WindowIndexConfig(cohort='TEST',output_root=self.base,chroms=['chr1'],samples=['S1'],
            sequence_variant_policy='snv_indel',window_size=16,offset_unit=1,blacklist_bed=bed))
        self.bam=self.base/'source.bam'
        bam_file(self.bam,self.sequence,[(self.sequence,[(0,n)],0,60,40),(self.hap2,[(0,6),(1,2),(0,14),(2,2),(0,n-22)],16,60,30)])
        self.directory=self.root/'compact'/LAYOUT/'chr1'

    def tearDown(self):self.temp.cleanup()

    def run_export(self):
        with patch.dict(os.environ,{'MFM_SAMPLE_CPU_BUDGET':'1','MFM_RNA_OUTPUT_MODE':'coverage_only'}):
            return export(self.root,self.bam,self.ref,self.ref,contigs=['chr1'],read_type='short_read',jobs=1)

    def test_export_preserves_windows_insertion_sequence_trace_reader_resume_cleanup(self):
        before=(self.root/'windows/windows.tsv').read_bytes();self.run_export()
        self.assertEqual(before,(self.root/'windows/windows.tsv').read_bytes())
        self.assertFalse(any(self.directory.glob('readset*')));self.assertFalse((self.directory/'depth.zarr').exists())
        self.assertEqual(validate_contig(self.directory,source=True)['HAP1']['status'],'pass')
        metrics={}
        with patch.dict(os.environ,{'MFM_PILEUP_VALIDATION_BACKEND':'aligned_pairs'}):
            self.assertEqual(validate_contig(self.directory,source=True,metrics=metrics)['HAP1']['status'],'pass')
        self.assertEqual(metrics['source_backend'],'aligned_pairs')
        group=zarr.open_group(str(self.directory/'pileup.zarr'),mode='r')
        from integrate_wgs_rna.pileup_features.legacy_export import ARRAYS
        expected_bytes=sum(group[hap][name].size*group[hap][name].dtype.itemsize for hap in ('HAP1','HAP2') for name in ARRAYS)
        self.assertEqual(metrics['pileup_array_bytes_read'],expected_bytes)
        q=self.root/'qc/final_package_validation.json';q.write_text(json.dumps(dict(sample='S1',status='pass',failed_checks=[])))
        reader=WindowReader(self.root);values=reader.load('chr1','HAP2',0)
        np.testing.assert_array_equal(values['dna_codes'][4:6],[0,3])
        np.testing.assert_array_equal(np.flatnonzero(values['pileup_insertion_mask']),[4,5])
        np.testing.assert_array_equal(values['wgs_depth'][4:6],[0,0])
        np.testing.assert_array_equal(values['pileup_features'][4:6,9],[.5,.5])
        np.testing.assert_array_equal(values['pileup_features'][4:6,4:9],np.zeros((2,5)))
        self.assertEqual(values['wgs_depth_coordinate_system'],'reference_projected_to_hap');reader.close()
        chunk=self.directory/'pileup.zarr/HAP1/features/0.0';mtime=chunk.stat().st_mtime_ns
        self.run_export();self.assertEqual(mtime,chunk.stat().st_mtime_ns)
        self._finalize()
        from handoff.cleanup import cleanup
        cleanup(self.root,'S1','wgs',output_root=self.base)
        self.assertFalse((self.root/'individual_fasta').exists())
        self.assertEqual(validate_contig(self.directory)['HAP2']['status'],'pass')
        self.assertEqual(WindowReader(self.root).load('chr1','HAP2',0)['pileup_features'].shape,(16,12))

    def _finalize(self):
        (self.root/'qc/wgs_preflight.json').write_text(json.dumps(dict(schema='mfm_preflight_qc_v2',sample='S1',status='pass')))
        events=[dict(schema='mfm_execution_provenance_v1',sample_id='S1',branch='WGS',stage=f'{i:02d}',event='skip',status='skipped',reason='실제 BAM/FASTA fixture 준비됨') for i in range(1,13)]
        (self.root/'metadata/execution_provenance.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in events))
        cmd=[sys.executable,str(ROOT/'scripts/validate_sample_package.py'),'--sample','S1','--root',str(self.root),'--mode','wgs',
             '--layout',self.directory.parent.name,'--chroms','chr1','--hap-read-mode','all_to_both','--wgs-payload','pileup_feature',
             '--output',str(self.root/'qc/final_package_validation.json'),'--completeness-tsv',str(self.root/'metadata/completeness.tsv')]
        proc=subprocess.run(cmd,capture_output=True,text=True)
        self.assertEqual(proc.returncode,0,proc.stdout+proc.stderr)

    def test_tampered_features_and_source_are_rejected(self):
        self.run_export();g=zarr.open_group(str(self.directory/'pileup.zarr'),mode='a');g['HAP1/features'][0,9]=.25
        with self.assertRaisesRegex(ValueError,'SHA'):validate_contig(self.directory)
        with self.assertRaises(ValueError):self.run_export()

    def test_native_source_is_independent_and_rejects_self_consistent_wrong_content(self):
        self.run_export()
        with patch.dict(os.environ,{'MFM_PILEUP_VALIDATION_BACKEND':'native'}):
            with patch('integrate_wgs_rna.pileup_features.compute.calculate',side_effect=RuntimeError('생산기를 호출하면 실패함')):
                validate_contig(self.directory,source=True)
            group=zarr.open_group(str(self.directory/'pileup.zarr'),mode='a')
            array=group['HAP1/features'];array[0,9]=.25
            path=self.directory/'contig_metadata.json';meta=json.loads(path.read_text())
            meta['pileup']['content_sha256']['HAP1']['features']=hashlib.sha256(array[:].tobytes()).hexdigest()
            path.write_text(json.dumps(meta))
            validate_contig(self.directory,source=False)
            with self.assertRaisesRegex(ValueError,'독립 reference source'):validate_contig(self.directory,source=True)

    def test_input_change_rejects_cached_pileup(self):
        self.run_export();stat=self.bam.stat()
        os.utime(self.bam,ns=(stat.st_atime_ns,stat.st_mtime_ns+1_000_000_000))
        with self.assertRaisesRegex(ValueError,'지문'):self.run_export()

    def test_interrupted_publication_is_resumable(self):
        rename=os.rename
        def fail_publish(source,destination):
            if Path(destination).name=='pileup.zarr':raise RuntimeError('시험용 공개 직전 중단')
            return rename(source,destination)
        with patch('integrate_wgs_rna.pileup_features.legacy_export.os.rename',side_effect=fail_publish):
            with self.assertRaisesRegex(RuntimeError,'중단'):self.run_export()
        self.assertFalse((self.directory/'pileup.zarr').exists())
        self.run_export()
        self.assertEqual(validate_contig(self.directory,source=True)['HAP1']['status'],'pass')

    def test_reference_and_chain_identity_fail_closed(self):
        wrong=self.base/'wrong.fa';wrong.write_text('>chr1\n'+'A'*len(self.sequence)+'\n');pysam.faidx(str(wrong))
        with self.assertRaisesRegex(Exception,'reference|M5'):
            export(self.root,self.bam,wrong,wrong,contigs=['chr1'],read_type='short_read',jobs=1)
        qc=self.root/'individual_fasta/personal_fasta_qc.json';data=json.loads(qc.read_text());data['reference_fasta']=str(wrong);qc.write_text(json.dumps(data))
        (qc.parent/MANIFEST).write_text(json.dumps(build(qc.parent,['chr1'],wrong)))
        with self.assertRaisesRegex(Exception,'reference|M5|무결성'):
            self.run_export()

    def test_pileup_reader_keeps_rna_window_signal(self):
        self.run_export()
        from integrate_wgs_rna.window_contract import read_windows
        rows=read_windows(self.directory/'windows.tsv',compact=True)
        for hap in ('HAP1','HAP2'):
            chosen=[r for r in rows if r['haplotype']==hap]
            lengths=np.array([int(r['actual_length']) for r in chosen])
            offsets=np.r_[0,np.cumsum(lengths)]
            g=zarr.open_group(str(self.root/f'rna_personal/panel/{hap}/coverage_windows.zarr'),mode='w')
            g.attrs.update(sample='S1',hap=hap)
            for name,value in dict(contig=np.array(['chr1']*len(chosen)),personal_start=np.array([int(r['raw_start']) for r in chosen]),
                personal_end=np.array([int(r['raw_start'])+int(r['actual_length']) for r in chosen]),offsets=offsets,
                depth=np.arange(int(offsets[-1]),dtype=np.uint32)).items():g.array(name,value)
        q=self.root/'qc/final_package_validation.json';q.write_text(json.dumps(dict(sample='S1',status='pass',failed_checks=[])))
        reader=WindowReader(self.root)
        try:
            value=reader.load('chr1','HAP2',0)
            np.testing.assert_array_equal(value['rna_depth'],np.arange(16,dtype=np.uint32))
            self.assertEqual(value['pileup_features'].shape,(16,12))
        finally:reader.close()

    def test_official_wgs_shell_skips_hap_alignment_and_finishes_cleanup(self):
        import shutil
        from handoff import planner
        target=self.root/'haplotagged/S1.illumina.hg38.haplotagged.bam';target.parent.mkdir()
        shutil.copy2(self.bam,target);shutil.copy2(Path(str(self.bam)+'.bai'),Path(str(target)+'.bai'))
        row=dict(sample_id='S1',platform='illumina',read_type='short_read',dataset_role='pre_training',project_id='TEST',cohort='TEST',
                 alignment_path=str(self.bam),raw_root=str(self.base),rna_path='',final_root=str(self.base),
                 work_root=str(self.base/'_work/test/S1'),common_vcf_dir='',persample_vcf_dir='',rna_bigwigs='',data_group='wgs',
                 source_decode_reference_fasta=str(self.ref),analysis_reference_fasta=str(self.ref),vcf_reference_fasta=str(self.ref),
                 vcf_source='illumina_panel',phase_source='provided')
        variants=self.base/'persample_vcf';variants.mkdir()
        row['persample_vcf_dir']=str(variants)
        vcf=variants/'S1.seq.vcf'
        pos=next(i for i,b in enumerate(self.sequence) if b in 'ACGT')
        ref=self.sequence[pos];alt='C' if ref=='A' else 'A'
        vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length='+str(len(self.sequence))+'>\n'
                       '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                       '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
                       'chr1\t'+str(pos+1)+'\t.\t'+ref+'\t'+alt+'\t60\tPASS\t.\tGT\t0|0\n')
        pysam.tabix_compress(str(vcf),str(vcf)+'.gz',force=True)
        pysam.tabix_index(str(vcf)+'.gz',preset='vcf',force=True)
        manifest=self.base/'manifest.tsv';planner.write_tsv(manifest,[row])
        (self.root/'metadata').mkdir(exist_ok=True);(self.root/'metadata/download_inputs.json').write_text('[]')
        (self.root/'qc').mkdir(exist_ok=True)
        events=[dict(schema='mfm_execution_provenance_v1',sample_id='S1',branch='WGS',stage=f'{i:02d}',event='skip',status='skipped',reason='시험용 준비 입력') for i in range(1,10)]
        (self.root/'metadata/execution_provenance.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in events))
        input_root=self.base/'input';input_root.mkdir();(input_root/'.mfm_handoff_store.json').write_text(json.dumps(dict(schema='mfm_handoff_store_v1')))
        (self.root/'qc/wgs_preflight.json').write_text(json.dumps(dict(schema='mfm_preflight_qc_v2',sample='S1',status='pass')))
        env={**os.environ,'PIPELINE_MANIFEST':str(manifest),'INPUT_ROOT':str(input_root),'OUTPUT_ROOT':str(self.base),'OUTPUT_ROOT_OVERRIDE':'0',
             'RUN_ID':'pileup_test','WGS_START_AT':'10','MFM_REPO_ROOT':str(ROOT/'vendor/MFM-DNARNAFM'),
             'MFM_CORE_ENV':str(Path(sys.executable).parent.parent),'PYTHONPATH':os.pathsep.join([str(ROOT),str(ROOT/'src'),str(ROOT/'vendor/MFM-DNARNAFM')]),
             'REF_SEQUENCE':str(self.ref),'REF_DECODE':str(self.ref),'ILLUMINA_REF_DECODE':str(self.ref),
             'MFM_SHORT_READ_OUTPUT_MODE':'pileup_feature','MFM_SEQUENCE_DEFER_RNA':'0','MFM_RNA_OUTPUT_MODE':'coverage_only',
             'CHROMS':'chr1','WINDOW_SIZE':'16','FIXED_WINDOW_START':'1','BLACKLIST_BED':str(self.base/'blacklist.bed'),
             'MFM_CPU_BUDGET':'2','MFM_SAMPLE_CPU_BUDGET':'2','MFM_CPU_POOL_DIR':str(self.base/'pool'),'MFM_CPU_POOL_STAGE_MAX':'2',
             'MFM_RAM_LIMIT_BYTES':str(128*1024**3),'CLEAN_INTERMEDIATES_ON_SUCCESS':'1','MFM_MODEL_ONLY_CLEANUP':'1',
             'MFM_EXTERNAL_CLEANUP_UNOWNED':'0','MFM_IO_BACKEND':'cpu','MFM_WGS_BACKEND':'cpu','MFM_RNA_BACKEND':'cpu',
             'MFM_BAM_BACKEND':'cpu','MFM_ZARR_BACKEND':'cpu','MFM_QC_BACKEND':'cpu','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
        env.pop('MFM_CPU_POOL_LEASE_FILE',None);env.pop('MFM_LONGREAD_MEMORY_PLAN',None)
        proc=subprocess.run(['bash','run_wgs_sample.sh','S1'],cwd=ROOT,env=env,capture_output=True,text=True,timeout=120)
        self.assertEqual(proc.returncode,0,proc.stdout[-7000:]+proc.stderr[-5000:])
        self.assertFalse((self.root/'realigned').exists());self.assertFalse((self.root/'individual_fasta').exists())
        self.assertEqual(json.loads((self.root/'qc/delivery_final_validation.json').read_text())['status'],'pass')
        events=[json.loads(line) for line in (self.root/'metadata/execution_provenance.jsonl').read_text().splitlines()]
        self.assertTrue(any(e.get('stage')=='10' and e.get('event')=='skip' for e in events))


if __name__=='__main__':unittest.main()
