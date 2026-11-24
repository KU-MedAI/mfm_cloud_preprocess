"""다운로드/WGS 의존성·reader의 없는 depth 의미·실제 내장 WGS 단계를 검증함."""
import csv,json,os,subprocess,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import pysam,zarr,numpy as np
from handoff import cli,planner
from handoff.resume import sample_downloads
from handoff.sequence_mode import selected,needs_wgs,validate_selection
from integrate_wgs_rna.window_reader import WindowReader

ROOT=Path(__file__).resolve().parents[1]

class SequenceModeTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name)
        self.config={'data_root':str(self.base/'input'),'output_root':str(self.base/'output'),'state_root':str(self.base/'state'),
                     'run_id':'sequence_test','pipeline_env':{'MFM_SHORT_READ_OUTPUT_MODE':'sequence_only'}}
        rows=planner.read_tsv(ROOT/'manifests/리드별/medai_short_read.tsv')
        self.row=next(r for r in rows if r['data_group']=='wgs');self.pair=next(r for r in rows if r['data_group']=='wgs_rna_pair')

    def tearDown(self):self.temp.cleanup()

    def test_wgs_only_raw_is_excluded_and_rna_evidence_is_explicit(self):
        resources=[{'kind':k,'path':k,'url':'fixture'} for k in ['wgs_alignment','wgs_index','illumina_panel','rna_alignment']]
        self.assertTrue(selected(self.config,self.row));self.assertFalse(needs_wgs(self.config,self.row))
        validate_selection(self.config,[self.row])
        self.assertEqual([r['kind'] for r in sample_downloads(self.config,self.row,resources)],['illumina_panel','rna_alignment'])
        self.assertTrue(needs_wgs(self.config,self.pair))
        with self.assertRaisesRegex(ValueError,'WGS 임시처리'):validate_selection(self.config,[self.pair])
        self.config['pipeline_env']['MFM_SEQUENCE_RNA_WGS_EVIDENCE']='temporary';validate_selection(self.config,[self.pair])
        self.assertEqual(resources,sample_downloads(self.config,self.pair,resources))
        long=dict(self.row,read_type='long_read',platform='hifi');self.assertFalse(selected(self.config,long));self.assertTrue(needs_wgs(self.config,long))

    def test_deferred_pair_filters_all_read_downloads_and_is_not_full_rna_completion(self):
        from handoff.sequence_mode import defer_rna,completion_mode,payload_matches,runtime_environment
        from handoff.rna_integrity import prepare_rna
        self.config['pipeline_env']['MFM_SEQUENCE_DEFER_RNA']='1'
        resources=[{'kind':k,'path':k,'url':'fixture'} for k in ['wgs_alignment','wgs_index','wgs_read','rna_alignment','rna_index','illumina_panel']]
        validate_selection(self.config,[self.pair]);self.assertTrue(defer_rna(self.config,self.pair))
        self.assertFalse(needs_wgs(self.config,self.pair));self.assertEqual(completion_mode(self.config,self.pair),'wgs')
        self.assertEqual([r['kind'] for r in sample_downloads(self.config,self.pair,resources)],['illumina_panel'])
        with patch('handoff.rna_integrity.fetch',side_effect=AssertionError('RNA를 수신하면 안 됨')):
            self.assertEqual(prepare_rna(self.config,self.pair,[]),(self.config,[]))
        self.assertEqual(runtime_environment(self.config,self.pair)['MFM_SEQUENCE_DEFER_RNA'],'1')
        root=self.base/'result';(root/'qc').mkdir(parents=True);p=root/'qc/final_package_validation.json'
        p.write_text(json.dumps({'wgs_payload':'sequence_only'}));self.assertFalse(payload_matches(self.config,self.pair,root))
        p.write_text(json.dumps({'wgs_payload':'sequence_only','completion_scope':'dna_sequence_only_rna_deferred'}))
        self.assertTrue(payload_matches(self.config,self.pair,root))
        full={**self.config,'pipeline_env':{'MFM_SHORT_READ_OUTPUT_MODE':'full'}}
        self.assertEqual(completion_mode(full,self.pair),'rna');self.assertTrue(needs_wgs(full,self.pair))
        self.config['pipeline_env']['MFM_SEQUENCE_DEFER_RNA']='bad'
        with self.assertRaises(ValueError):validate_selection(self.config,[self.pair])

    def test_bundled_only_admission_keeps_ram_guard(self):
        from handoff.sequence_mode import uses_bundled_only
        from handoff.queue_resources import download_admit
        assets=self.base/'assets';assets.mkdir();(assets/'panel.gz').write_bytes(b'VCF')
        (assets/'manifest.json').write_text(json.dumps({'schema':'mfm_bundled_assets_v1','files':[{'path':'panel.gz','size_bytes':3}]}))
        self.config['bundled_assets_root']=str(assets)
        resources=[{'kind':'illumina_panel','shared':True,'path':'panel.gz'}]
        with patch('handoff.downloads.sample_resources',return_value=resources):
            self.assertTrue(uses_bundled_only(self.config,self.row,ROOT))
            self.assertFalse(uses_bundled_only(self.config,self.pair,ROOT))
            (assets/'panel.gz').write_bytes(b'damaged');self.assertFalse(uses_bundled_only(self.config,self.row,ROOT))
        pool=self.base/'pool';pool.mkdir();(pool/'policy.json').write_text(json.dumps({'max_download_jobs':1,'download_min_free_bytes':10**18}))
        leases={'existing':[{'branch':'DOWNLOAD','memory_bytes':1_000_000_000}]}
        with patch('handoff.queue_resources.live_leases',return_value=leases),patch('handoff.queue_resources.cgroup_limits',return_value=(0,None,None)),patch('handoff.queue_resources.effective_memory',return_value=10**12):
            self.assertFalse(download_admit(self.config,pool)[0])
            self.assertTrue(download_admit(self.config,pool,local_only=True)[0])
            with patch('handoff.queue_resources.effective_memory',return_value=1_000_000_000):
                self.assertFalse(download_admit(self.config,pool,local_only=True)[0])

    def test_vcf_workers_use_granted_cpu_and_spread_across_samples(self):
        from integrate_wgs_rna.vcf_execution import worker_affinity
        with patch('os.sched_getaffinity',return_value=set(range(24,96))):
            with patch.dict(os.environ,{'MFM_WGS_PAYLOAD':'full'}):
                self.assertEqual(worker_affinity('A',2,0),[24,25])
                self.assertEqual(worker_affinity('B',2,1),[26,27])
            with patch.dict(os.environ,{'MFM_WGS_PAYLOAD':'sequence_only'}):
                grants=[worker_affinity('S'+str(i),2,0) for i in range(20)]
                self.assertTrue(all(len(g)==2 and len(set(g))==2 and set(g)<=set(range(24,96)) for g in grants))
                self.assertGreater(len({c for g in grants for c in g}),20)
                self.assertEqual(grants[0],worker_affinity('S0',2,0))

    def test_completed_full_payload_is_not_reused_as_sequence(self):
        from types import SimpleNamespace
        from handoff.local_queue import Adapter
        from handoff.sequence_mode import payload_matches
        root=self.base/'S';(root/'qc').mkdir(parents=True);qc=root/'qc/final_package_validation.json'
        adapter=Adapter.__new__(Adapter);adapter.config=self.config
        adapter.api=SimpleNamespace(sample_root=lambda *a:root,processing_complete=lambda *a:True)
        qc.write_text(json.dumps({'status':'pass'}))
        self.assertFalse(payload_matches(self.config,self.row,root));self.assertFalse(adapter.verified(self.row))
        qc.write_text(json.dumps({'status':'pass','wgs_payload':'sequence_only'}))
        self.assertTrue(adapter.verified(self.row))
        full={**self.config,'pipeline_env':{'MFM_SHORT_READ_OUTPUT_MODE':'full'}}
        self.assertFalse(payload_matches(full,self.row,root))

    def test_binding_without_wgs_does_not_read_nonexistent_cram(self):
        refs={k:str(self.base/(k+'.fa')) for k in ['autosomes','ncbi','broad','illumina','gtf','rna_extra']}
        with patch('handoff.references.exact_cram_reference',side_effect=AssertionError('WGS를 읽으면 안 됨')),patch('handoff.references.variant_reference',return_value=refs['autosomes']),patch('handoff.bundled_assets.panel_directory',return_value=self.base/'panel'):
            manifest,bound=cli.bind_manifest(self.config,self.row,[],refs)
        self.assertEqual(bound['wgs_payload'],'sequence_only');self.assertTrue(manifest.is_file())

    def test_actual_wgs_shell_preserves_window_policy_and_cleans_without_raw(self):
        self.run_actual_sequence(False)

    def test_actual_pair_shell_preserves_identity_defers_rna_and_keeps_unused_raw(self):
        self.run_actual_sequence(True)

    def test_actual_sequence_with_fixed_cpu_keeps_same_window_and_qc_contract(self):
        self.run_actual_sequence(False,fixed=True)

    def test_actual_public_long_sequence_without_wgs_reads(self):
        self.run_actual_sequence(False,long=True,preserve_raw=True)

    def run_actual_sequence(self,deferred,fixed=False,long=False,preserve_raw=False):
        from integrate_wgs_rna.sequence_only import LAYOUT
        base=self.base;sequence='NNN'+'ACGT'*25+'A'*45+'N'+'A'*54+'N'*100+'A'*12+'N'+'A'*12+'NNNN'
        ref=base/'reference.fa';ref.write_text('>chr1\n'+sequence+'\n');pysam.faidx(str(ref))
        panel=base/'panel';panel.mkdir();vcf=panel/'TEST.chr1.vcf'
        # chr1:4의A→C 이형접합 변이로 실제03/09를 실행함.
        vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length='+str(len(sequence))+'>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t4\t.\tA\tC\t60\tPASS\t.\tGT\t0|1\nchr1\t8\t.\tA\tAT\t60\tPASS\t.\tGT\t1|1\n')
        pysam.tabix_compress(str(vcf),str(vcf)+'.gz',force=True);pysam.tabix_index(str(vcf)+'.gz',preset='vcf',force=True)
        bed=base/'blacklist.bed';bed.write_text('chr1\t3\t103\n')
        row=dict(self.row,sample_id='S1',cohort='TEST',project_id='TEST',source_sample_id='S1',alignment_path=str(base/'NEVER_DOWNLOAD.cram'),
                 common_vcf_dir=str(panel),persample_vcf_dir=str(base/'output/_work/sequence_test/S1/persample_vcf'),
                 source_decode_reference_fasta=str(ref),analysis_reference_fasta=str(ref),vcf_reference_fasta=str(ref),rna_path='',
                 final_root='',work_root='',raw_root='',public_vcf_path='',public_vcf_sample_id='')
        if long:row.update(alignment_path='',read_type='long_read',platform='hifi',vcf_source='public',phase_source='provided',public_vcf_path=str(vcf)+'.gz',public_vcf_sample_id='S1',vcf_callset='PAV_v2.4.0.1')
        if deferred:row.update(data_group='wgs_rna_pair',rna_path=str(base/'NEVER_DOWNLOAD_RNA.bam'))
        manifest=base/'manifest.tsv';planner.write_tsv(manifest,[row])
        sample=base/'output'/row['data_group']/'TEST/S1';(sample/'metadata').mkdir(parents=True);(sample/'metadata/download_inputs.json').write_text('[]')
        data=base/'input';data.mkdir();(data/'.mfm_handoff_store.json').write_text(json.dumps({'schema':'mfm_handoff_store_v1'}))
        env={**os.environ,'PIPELINE_MANIFEST':str(manifest),'INPUT_ROOT':str(base/'input'),'OUTPUT_ROOT':str(base/'output'),
             'OUTPUT_ROOT_OVERRIDE':'1','RUN_ID':'sequence_test','MFM_REPO_ROOT':str(ROOT/'vendor/MFM-DNARNAFM'),
             'MFM_CORE_ENV':str(Path(os.sys.executable).parent.parent),'PYTHONPATH':os.pathsep.join([str(ROOT),str(ROOT/'src'),str(ROOT/'vendor/MFM-DNARNAFM')]),
             'REF_SEQUENCE':str(ref),'ILLUMINA_REF_DECODE':str(ref),'REF_DECODE':str(ref),'PANEL_BASE':str(panel),'PANEL_LOCAL_DIR':str(panel),
             'CHROMS':'chr1','WINDOW_SIZE':'100','FIXED_WINDOW_START':'1','BLACKLIST_BED':str(bed),
             'MFM_SHORT_READ_OUTPUT_MODE':'sequence_only','MFM_CPU_BUDGET':'2','MFM_SAMPLE_CPU_BUDGET':'2',
             'MFM_CPU_POOL_DIR':str(base/'pool'),'MFM_CPU_POOL_STAGE_MAX':'2','MFM_RAM_LIMIT_BYTES':str(128*1024**3),
             'CLEAN_INTERMEDIATES_ON_SUCCESS':'1','MFM_MODEL_ONLY_CLEANUP':'1','MFM_EXTERNAL_CLEANUP_UNOWNED':'0',
             'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MFM_ZARR_BACKEND':'cpu','MFM_QC_BACKEND':'cpu',
             'MFM_WGS_BACKEND':'cpu','MFM_RNA_BACKEND':'cpu','MFM_BAM_BACKEND':'cpu','MFM_IO_BACKEND':'cpu'}
        if preserve_raw:env['MFM_MODEL_ONLY_CLEANUP']='0'
        if long:env['MFM_LONG_READ_OUTPUT_MODE']='sequence_only'
        if fixed:env['MFM_CPU_ALLOCATION_MODE']='fixed_sample'
        env['MFM_SEQUENCE_DEFER_RNA']='1' if deferred else '0'
        preserved=[]
        if deferred or preserve_raw:
            import hashlib
            for area in ['wgs','rna']:
                f=data/area/'TEST/S1/original.bam';f.parent.mkdir(parents=True,exist_ok=True);f.write_bytes(b'preserve unused raw')
                records=data/'.download_records';records.mkdir(exist_ok=True)
                key=hashlib.sha256(str(f.resolve()).encode()).hexdigest();(records/(key+'.json')).write_text(json.dumps({'owned':True}))
                preserved.append(f)
        env.pop('MFM_CPU_POOL_LEASE_FILE',None);env.pop('MFM_LONGREAD_MEMORY_PLAN',None)
        proc=subprocess.run(['bash','run_sample.sh' if deferred else 'run_wgs_sample.sh','S1'],cwd=ROOT,env=env,capture_output=True,text=True,timeout=90)
        self.assertEqual(proc.returncode,0,proc.stdout[-3000:]+proc.stderr[-5000:])
        self.assertFalse((base/'NEVER_DOWNLOAD.cram').exists());self.assertFalse((sample/'individual_fasta').exists())
        if deferred:
            meta=json.loads((sample/'metadata/sample.json').read_text());qc=json.loads((sample/'qc/final_package_validation.json').read_text())
            self.assertTrue(meta['has_rna']);self.assertEqual(meta['data_group'],'wgs_rna_pair')
            self.assertEqual(meta['rna_processing_status'],'deferred');self.assertEqual(qc['mode'],'wgs')
            self.assertEqual(qc['completion_scope'],'dna_sequence_only_rna_deferred')
            self.assertFalse((sample/'rna_personal').exists());self.assertFalse((base/'NEVER_DOWNLOAD_RNA.bam').exists())
            self.assertTrue(all(f.read_bytes()==b'preserve unused raw' for f in preserved))
            provenance=[json.loads(line) for line in (sample/'metadata/execution_provenance.jsonl').read_text().splitlines()]
            self.assertFalse(any(x['branch']=='RNA' for x in provenance))
            for stage in ['01','02','06','10' if (ROOT/'09_personal_fasta.sh').exists() else '08']:
                self.assertTrue(any(x['branch']=='WGS' and x['stage']==stage and x['event']=='skip' for x in provenance))
        if preserve_raw:self.assertTrue(all(f.read_bytes()==b'preserve unused raw' for f in preserved))
        self.assertEqual(json.loads((sample/'qc/final_package_validation.json').read_text())['status'],'pass')
        self.assertEqual(json.loads((sample/'qc/delivery_final_validation.json').read_text())['status'],'pass')
        self.assertFalse(any((sample/'compact'/LAYOUT).rglob('depth.zarr')))
        reader=WindowReader(sample);data=reader.load('chr1','HAP1',0)
        with (sample/'metadata/wgs_windows.tsv').open() as handle:
            retained=[r for r in csv.DictReader(handle,delimiter='\t') if r['haplotype']=='HAP1' and r['kept']=='1']
        self.assertIsNone(data['wgs_depth']);self.assertEqual(reader.window_count('chr1','HAP1'),len(retained))
        reader.close()

if __name__=='__main__':unittest.main()
