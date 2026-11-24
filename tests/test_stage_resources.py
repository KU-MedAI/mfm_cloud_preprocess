"""public HiFi 단계예측·legacy회귀·실제lease투입계약을검증함."""
from pathlib import Path
from unittest.mock import patch
import tempfile
import unittest
from handoff.stage_resources import longread_memory_plan,starting_memory
from handoff.queue_resources import full_ram,active_ram,GIB


class Tests(unittest.TestCase):
    def row(self,**kwargs):
        return dict(sample_id='S',platform='hifi',read_type='long_read',vcf_source='public',phase_source='provided',data_group='wgs',known_input_total_bytes='488000000000')|kwargs

    def test_입력바이트크기아닌reference와buffer로예약(self):
        c={'sample_cpus':68};a=longread_memory_plan(self.row(),c);b=longread_memory_plan(self.row(known_input_total_bytes='1000000000'),c)
        self.assertEqual(a['admission_bytes'],b['admission_bytes']);self.assertLess(a['admission_bytes'],900_000_000_000)
        self.assertGreaterEqual(a['admission_bytes'],max(a['stage_estimates_bytes'].values()))

    def test_assembly나selfcall미검증경로는보수식유지(self):
        self.assertIsNone(longread_memory_plan(self.row(vcf_source='clair3')))
        self.assertIsNone(longread_memory_plan(self.row(),{'pipeline_env':{'SEQUENCE_SOURCE':'denovo'}}))
        legacy=full_ram([self.row()],{'pipeline_env':{'MFM_LONGREAD_MEMORY_MODEL':'legacy'}})
        self.assertGreater(legacy,900_000_000_000)

    def test_RNA별도예산과작은CPU메모리정책(self):
        plain=longread_memory_plan(self.row());rna=longread_memory_plan(self.row(data_group='wgs_rna_pair'))
        self.assertGreaterEqual(rna['admission_bytes'],max(rna['stage_estimates_bytes'].values()))
        self.assertIn('rna_star_and_sort',rna['stage_estimates_bytes'])

    def test_단계예약87GB를100GiB로부풀리지않음(self):
        row=self.row(platform='illumina',read_type='short_read')
        with patch('handoff.queue_resources.live_leases',return_value={'S':[{'branch':'WGS','stage':'10','memory_bytes':87_000_000_000}]}):
            self.assertEqual(active_ram(Path('/tmp'),[row]),16*GIB+87_000_000_000)

    def test_lease없거나예약0이면보수모델복귀(self):
        row=self.row()
        for leases in [{},{'S':[{'branch':'WGS','memory_bytes':0}]}]:
            with patch('handoff.queue_resources.live_leases',return_value=leases):
                self.assertEqual(active_ram(Path('/tmp'),[row]),16*GIB+starting_memory(row))

    def test_후보자자신은기존lease가있어도전체시작예산(self):
        row=self.row()
        with patch('handoff.queue_resources.live_leases',return_value={'S':[{'branch':'WGS','memory_bytes':4*GIB}]}):
            self.assertEqual(active_ram(Path('/tmp'),[],row),16*GIB+starting_memory(row))


    def test_실제단계예약과모델예측이양쪽layout에서동일(self):
        import importlib.util,json
        spec=importlib.util.spec_from_file_location('hifi_pool',Path(__file__).resolve().parents[1]/'scripts/cpu_pool.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        plan=longread_memory_plan(self.row(data_group='wgs_rna_pair'),{'sample_cpus':68})
        self.assertGreaterEqual(plan['stage_estimates_bytes']['consensus'],4*GIB*35)
        for layout,stages in [('original',{'01':'prepare_alignment','09':'consensus_entry','10':'personal_realign_entry','12':'compact_entry'}),('legacy',{'01':'prepare_alignment','07':'consensus_entry','08':'personal_realign_entry','10':'compact_entry'})]:
            with patch.dict('os.environ',{'MFM_LONGREAD_MEMORY_PLAN':json.dumps(plan),'MFM_CPU_POOL_LAYOUT':layout}):
                for stage,label in stages.items():self.assertEqual(module.stage_memory('WGS',stage),plan['stage_estimates_bytes'][label])
                self.assertEqual(module.stage_memory('RNA','13' if layout=='legacy' else '15'),plan['stage_estimates_bytes']['rna_star_index'])
                self.assertEqual(module.stage_memory('RNA','14' if layout=='legacy' else '16'),plan['stage_estimates_bytes']['rna_star_align'])


    def test_VCF입장과haplotag예산분리및구형계획호환(self):
        import importlib.util,json
        from integrate_wgs_rna.public_vcf_parallel import entry_memory
        spec=importlib.util.spec_from_file_location('vcf_memory_pool',Path(__file__).resolve().parents[1]/'scripts/cpu_pool.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        plan=longread_memory_plan(self.row(),{'sample_cpus':68})
        self.assertEqual(plan['stage_estimates_bytes']['public_vcf_prepare'],entry_memory(68))
        self.assertEqual(plan['stage_estimates_bytes']['haplotag'],128*GIB)
        for layout in ['original','legacy']:
            for older in [False,True]:
                candidate=json.loads(json.dumps(plan))
                if older:
                    candidate['stage_estimates_bytes'].pop('public_vcf_prepare')
                    candidate['stage_estimates_bytes'].pop('haplotag')
                    candidate['stage_estimates_bytes']['public_vcf_and_haplotag']=128*GIB
                with patch.dict('os.environ',{'MFM_LONGREAD_MEMORY_PLAN':json.dumps(candidate),'MFM_CPU_POOL_LAYOUT':layout,'MFM_CPU_BUDGET':'68'}):
                    self.assertEqual(module.stage_memory('WGS','03'),entry_memory(68))
                    self.assertEqual(module.stage_memory('WGS','06'),128*GIB)
        self.assertEqual(starting_memory(self.row(),{'sample_cpus':68,'pipeline_env':{'WGS_START_AT':'03'}}),entry_memory(68))
        self.assertEqual(starting_memory(self.row(),{'pipeline_env':{'WGS_START_AT':'06'}}),128*GIB)

    def test_공유풀실제실행에서구형128GiB계획으로작은RAM입장(self):
        import json,os,subprocess,sys
        from integrate_wgs_rna.public_vcf_parallel import entry_memory
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);plan=longread_memory_plan(self.row(),{'sample_cpus':68})
            plan['stage_estimates_bytes'].pop('public_vcf_prepare');plan['stage_estimates_bytes'].pop('haplotag')
            plan['stage_estimates_bytes']['public_vcf_and_haplotag']=128*GIB
            result_path=root/'child.json'
            code="import os,json,pathlib;pathlib.Path(os.environ['CHECK_RESULT']).write_text(pathlib.Path(os.environ['MFM_CPU_POOL_LEASE_FILE']).read_text())"
            env={**os.environ,'MFM_LONGREAD_MEMORY_PLAN':json.dumps(plan),'MFM_CPU_BUDGET':'8','MFM_RAM_LIMIT_BYTES':str(20*GIB),'MFM_CPU_POOL_LAYOUT':'original','CHECK_RESULT':str(result_path),'MFM_CPU_UTILIZATION':'maximum'}
            command=[sys.executable,str(Path(__file__).resolve().parents[1]/'scripts/cpu_pool.py'),'--pool-dir',str(root/'pool'),'--total-cpus','8','--branch','WGS','--stage','03','--sample','S','--',sys.executable,'-c',code]
            run=subprocess.run(command,env=env,capture_output=True,text=True,timeout=15)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            self.assertEqual(json.loads(result_path.read_text())['memory_bytes'],entry_memory(8))
            self.assertFalse(list((root/'pool/leases').glob('*.json')))

    def test_FASTA입장은과거계획과payload에도4GiB(self):
        import importlib.util,json,os
        spec=importlib.util.spec_from_file_location('fasta_entry_pool',Path(__file__).resolve().parents[1]/'scripts/cpu_pool.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        plan=longread_memory_plan(self.row(),{'sample_cpus':68})
        self.assertEqual(plan['stage_estimates_bytes']['consensus_entry'],4*GIB)
        self.assertGreaterEqual(plan['stage_estimates_bytes']['consensus'],140*GIB)
        for layout,stage in [('original','09'),('legacy','07')]:
            for declared in [plan,{'schema':'mfm_longread_stage_memory_v1','stage_estimates_bytes':{'consensus':128*GIB}},None]:
                for payload in ['full','sequence_only']:
                    env={'MFM_CPU_POOL_LAYOUT':layout,'MFM_WGS_PAYLOAD':payload,'MFM_LONGREAD_MEMORY_PLAN':json.dumps(declared) if declared else ''}
                    with patch.dict(os.environ,env):self.assertEqual(module.stage_memory('WGS',stage),4*GIB)
        self.assertEqual(starting_memory(self.row(),{'pipeline_env':{'WGS_START_AT':'09'}}),4*GIB)

    def test_FASTA계산은8GiB없으면시작하지않음(self):
        import json,os
        from integrate_wgs_rna.cpu_executor import WorkerBudget
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();(root/'waiters').mkdir()
            (root/'policy.json').write_text(json.dumps({'live_limits':True,'cpu_pool_limit':8,'memory_limit_bytes':20*GIB}))
            own=root/'leases/own.json';own.write_text(json.dumps({'cpus':2,'total_cpus':8,'memory_bytes':4*GIB}))
            other=root/'leases/other.json';other.write_text(json.dumps({'cpus':6,'memory_bytes':14*GIB}))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(own),'MFM_RAM_LIMIT_BYTES':str(20*GIB),'MFM_MEMORY_STAGE_BASE_BYTES':str(4*GIB)}):
                budget=WorkerBudget(1,44,True,4*GIB,cpus_per_worker=2)
                self.assertEqual(budget.capacity(44),0)
                self.assertEqual(json.loads(own.read_text())['memory_bytes'],4*GIB)
                other.write_text(json.dumps({'cpus':6,'memory_bytes':10*GIB}))
                self.assertEqual(budget.capacity(44),1)
                self.assertEqual(json.loads(own.read_text())['memory_bytes'],8*GIB)

    def test_과거128GiB계획도실제FASTA입장은4GiB(self):
        import json,os,subprocess,sys
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);result=root/'child.json'
            plan={'schema':'mfm_longread_stage_memory_v1','stage_estimates_bytes':{'consensus':128*GIB}}
            env={**os.environ,'MFM_LONGREAD_MEMORY_PLAN':json.dumps(plan),'MFM_CPU_BUDGET':'8',
                 'MFM_RAM_LIMIT_BYTES':str(16*GIB),'MFM_CPU_POOL_LAYOUT':'original',
                 'MFM_WGS_PAYLOAD':'full','CHECK_RESULT':str(result),'MFM_CPU_UTILIZATION':'maximum'}
            code="import os,pathlib;pathlib.Path(os.environ['CHECK_RESULT']).write_text(pathlib.Path(os.environ['MFM_CPU_POOL_LEASE_FILE']).read_text())"
            command=[sys.executable,str(Path(__file__).resolve().parents[1]/'scripts/cpu_pool.py'),
                     '--pool-dir',str(root/'pool'),'--total-cpus','8','--branch','WGS','--stage','09',
                     '--sample','S','--',sys.executable,'-c',code]
            child=subprocess.run(command,env=env,capture_output=True,text=True,timeout=15)
            self.assertEqual(child.returncode,0,child.stdout+child.stderr)
            self.assertEqual(json.loads(result.read_text())['memory_bytes'],4*GIB)
            self.assertFalse(list((root/'pool/leases').glob('*.json')))


if __name__=='__main__':unittest.main()
