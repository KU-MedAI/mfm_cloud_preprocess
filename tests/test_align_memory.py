"""정렬 메모리 추정·HAP 경계 재배정·실행중 보정의 안전 계약을 검증함."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import json
import os
import tempfile
import unittest

from integrate_wgs_rna.memory_policy import GIB
from integrate_wgs_rna.realign.align_memory import alignment_memory
from integrate_wgs_rna.realign.align_resources import alignment_plan
from integrate_wgs_rna.realign.phases import RealignPhases,_atomic
from integrate_wgs_rna.realign.live_budget import propose,apply

AUTO={'MFM_WGS_ALIGN_REQUESTED':'auto','MFM_WGS_SORT_REQUESTED':'auto','MFM_WGS_HAP_REQUESTED':'auto'}


class Tests(unittest.TestCase):
    def test_1HAP예약이192GiB고정아니며_산식여유포함(self):
        m=alignment_memory(SimpleNamespace(sort_mem='2G'),{1:3_000_000_000},1,8,2,{})
        self.assertLess(m['reserved_memory_bytes'],192*GIB)
        self.assertGreater(m['reserved_memory_bytes'],3_000_000_000*16+24*GIB+6*GIB)
        self.assertEqual(m['reserved_memory_bytes']%GIB,0)

    def test_스레드_reference_sort_HAP증가를모두예약에반영(self):
        c=SimpleNamespace(sort_mem='2G');lengths={1:3_000_000_000,2:3_100_000_000}
        baseline=alignment_memory(c,lengths,1,4,0,{})['reserved_memory_bytes']
        for h,t,s,l in [(1,32,0,lengths),(1,4,4,lengths),(2,4,0,lengths),(1,4,0,{1:4_000_000_000})]:
            self.assertGreater(alignment_memory(c,l,h,t,s,{})['reserved_memory_bytes'],baseline)

    def test_samtools_sort는_0과1스레드에서_버퍼1개만예약(self):
        config = SimpleNamespace(sort_mem='64G', platform='hifi')
        for threads, buffers in [(0, 1), (1, 1), (4, 4)]:
            model = alignment_memory(config, {1:3_000_000_000}, 1, 4, threads, {})
            self.assertEqual(model['sort_buffers_bytes'], 64*GIB*buffers)

    def test_legacy모델선택과설정오류(self):
        c=SimpleNamespace(sort_mem='2G')
        self.assertEqual(alignment_memory(c,{1:1000},1,1,0,{'MFM_REALIGN_MEMORY_MODEL':'legacy192'})['reserved_memory_bytes'],192*GIB)
        with self.assertRaises(ValueError):alignment_memory(c,{1:1000},1,1,0,{'MFM_REALIGN_MEMORY_MODEL':'bad'})

    def test_RAM이192GiB미만이어도증명된reference계획으로시작(self):
        from test_align_resources import Plans
        p=alignment_plan(Plans().config(),{1:3_000_000_000},12,110*GIB,AUTO)
        self.assertIsNotNone(p);self.assertLessEqual(p['reserved_memory_bytes'],110*GIB)

    def test_실제BWA_HAP1다음반환된CPU로HAP2확장_QC동일(self):
        from test_stage22_bwa import Stage22BwaIntegrationTest,module
        original=module.run_personal_realign_generation;reserve=RealignPhases.align_config;seen=[];calls=[]
        def allocate(phase,config,lengths):
            calls.append(lengths);phase.total=12 if len(calls)<=2 else 24
            return reserve(phase,config,lengths)
        def observe(config):
            result=original(config)
            p=config.output_root/'TEST/qc/TEST.stage11_timing.json'
            if p.exists():seen.append(json.loads(p.read_text()))
            return result
        with patch.dict(os.environ,dict(AUTO,MFM_CPU_BUDGET='24')),patch.object(RealignPhases,'align_config',allocate),patch.object(module,'run_personal_realign_generation',side_effect=observe):
            Stage22BwaIntegrationTest('test_pair_singleton_duplicate_and_chain_contract').debug()
        self.assertEqual(seen[0]['effective_threads_by_hap']['HAP1']['align_per_hap'],7)
        self.assertEqual(seen[0]['effective_threads_by_hap']['HAP2']['align_per_hap'],17)

    def test_완료HAP의CPU_RAM반환_남은HAP은그대로(self):
        from test_align_resources import Plans
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'leases').mkdir();(root/'waiters').mkdir();owner=root/'leases/a.json'
            _atomic(owner,{'cpus':22,'total_cpus':68,'memory_bytes':192*GIB,'manager_pid':os.getpid(),'realign_phase':'align_restore_sort_qc','effective_alignment':{'hap_jobs':2}})
            config=Plans().config(root);config.threads=6;config.sort_threads=2
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(owner)},clear=True):
                RealignPhases().hap_completed(config,{2:3_000_000_000})
            row=json.loads(owner.read_text());self.assertEqual(row['cpus'],11);self.assertLess(row['memory_bytes'],192*GIB);self.assertEqual(row['effective_alignment']['hap_jobs'],1)

    def test_QC는완료HAP제외한실제동시수를사용(self):
        from integrate_wgs_rna.realign.output_scan import output_scan_plan
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);output=root/'a.bam';output.write_bytes(b'0');lease=root/'lease.json';lease.write_text(json.dumps({'cpus':8,'memory_bytes':80*GIB,'effective_alignment':{'hap_jobs':1}}))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(lease)},clear=True):
                p=output_scan_plan(SimpleNamespace(jobs=1,hap_jobs=2,sort_threads=1),output,[('chr1',1000)],True,None)
            self.assertEqual(p['cpu_budget'],8);self.assertEqual(p['concurrent_hap_samples'],1)


class LiveTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name);(self.root/'leases').mkdir();self.lease=self.root/'leases/a.json'
        self.row={'sample':'S','branch':'WGS','stage':'10','manager_pid':123,'cpus':6,'memory_bytes':192*GIB,'realign_phase':'align_restore_sort_qc','effective_alignment':{'hap_jobs':1,'threads':3,'sort_threads':0}}
        self.lease.write_text(json.dumps(self.row));self.cohort=self.root/'cohort';sample=self.cohort/'S'
        for hap in [1,2]:
            d=sample/'individual_fasta'/f'HAP{hap}';d.mkdir(parents=True)
            for n in range(1,23):(d/f'chr{n}.fa.fai').write_text(f'chr{n}\t137000000\t0\t60\t61\n')
        logs=sample/'logs/personal_realign';logs.mkdir(parents=True)
        (logs/'S.HAP1.realign.sh').write_text(' | "$SAMTOOLS" sort -@ 0 -m 2G -o "$PART" -\n')
        self.tree={'manager_start':456,'rss_bytes':20*GIB,'bwa':[{'pid':999,'start':789,'argv':['bwa-mem2','mem','-t','3']}]}
        self.patch=patch('integrate_wgs_rna.realign.live_budget.tree_snapshot',side_effect=lambda pid:self.tree);self.patch.start();self.addCleanup(self.patch.stop)

    def test_검증후메모리만감소_CPU_PID유지(self):
        p=propose(self.lease,self.cohort);self.assertLess(p['new_memory_bytes'],192*GIB)
        r=apply(p);self.assertEqual(r['status'],'applied');row=json.loads(self.lease.read_text());self.assertEqual(row['cpus'],6);self.assertEqual(row['manager_pid'],123)

    def hifi_fixture(self):
        self.row.update(cpus=22, memory_bytes=450*GIB,
            effective_alignment={'hap_jobs':2,'threads':4,'sort_threads':1,'input_decode_threads':1},
            hap_resources={str(h):{'cpus':11,'memory_bytes':225*GIB,'phase':'align_restore_sort_qc'} for h in [1,2]})
        self.lease.write_text(json.dumps(self.row))
        logs=self.cohort/'S/logs/personal_realign'
        for hap in [1,2]:
            (logs/f'S.HAP{hap}.realign.sh').write_text(' | "$SAMTOOLS" sort -@ 1 -m 64G -o "$PART" -\n')
        self.tree.update(rss_bytes=140*GIB, bwa=[],
            minimap2=[{'pid':900+h,'start':700+h,'argv':['minimap2','-ax','map-hifi','-t','4']} for h in [1,2]],
            sort=[{'pid':910+h,'start':710+h,'argv':['samtools','sort','-@','1','-m','64G']} for h in [1,2]])

    def test_HiFi도_PID와sort검증후_예약만줄이고_HAP별몫보존(self):
        self.hifi_fixture();proposal=propose(self.lease,self.cohort);apply(proposal)
        row=json.loads(self.lease.read_text())
        self.assertLess(row['memory_bytes'],450*GIB)
        self.assertGreaterEqual(row['memory_bytes'],280*GIB)
        self.assertEqual(row['cpus'],22)
        self.assertEqual(sum(s['memory_bytes'] for s in row['hap_resources'].values()),row['memory_bytes'])
        self.assertEqual([s['cpus'] for s in row['hap_resources'].values()],[11,11])
        self.assertEqual(row['effective_alignment']['memory_model']['sort_buffers_bytes'],128*GIB)

    def test_HiFi실제sort버퍼가다르면보정거부(self):
        self.hifi_fixture();self.tree['sort'][0]['argv'][-1]='32G'
        with self.assertRaises(ValueError):propose(self.lease,self.cohort)

    def test_HiFi도구_PID변경시적용거부(self):
        self.hifi_fixture();proposal=propose(self.lease,self.cohort)
        self.tree['minimap2'][0]['start']+=1
        with self.assertRaises(ValueError):apply(proposal)

    def test_lease가변경되면보정거부(self):
        p=propose(self.lease,self.cohort);self.row['cpus']=5;self.lease.write_text(json.dumps(self.row))
        with self.assertRaises(ValueError):apply(p)

    def test_PID가재사용되면보정거부(self):
        p=propose(self.lease,self.cohort);self.tree['manager_start']=999
        with self.assertRaises(ValueError):apply(p)

    def test_reference가변경되면보정거부(self):
        p=propose(self.lease,self.cohort);f=self.cohort/'S/individual_fasta/HAP1/chr1.fa.fai';f.write_text('chr1\t137000001\t0\t60\t61\n')
        with self.assertRaises(ValueError):apply(p)

    def test_RSS가높으면감소안함(self):
        self.tree['rss_bytes']=100*GIB
        with self.assertRaises(ValueError):propose(self.lease,self.cohort)

    def test_index나완료상태에는보정안함(self):
        self.row['realign_phase']='bwa_index';self.lease.write_text(json.dumps(self.row))
        with self.assertRaises(ValueError):propose(self.lease,self.cohort)


if __name__=='__main__':unittest.main()
