"""필드·순서·전수해시·자원 상한을 실제 파일로 대조함."""
import hashlib,json,os,subprocess,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import pysam
from integrate_wgs_rna.compact.trace_hash import source_digests,digests,segment_digests
from integrate_wgs_rna.gpu_cram import encode,plan
from integrate_wgs_rna.gpu_stage_resources import StageBudget
from integrate_wgs_rna.realign.hap_resources import HaplotypeResources
from integrate_wgs_rna.realign.phases import _atomic
from integrate_wgs_rna.memory_policy import GIB
from integrate_wgs_rna.realign.align_resources import alignment_plan


class BottleneckTests(unittest.TestCase):
    def test_hash_is_exact_old_json_and_detects_changes(self):
        values=[('ACGTN',1,9223372036854775807,-1,500,0,(16,32)),('RYSWKMBDHVN',0,0,100,999,1,()),('',2,42,0,0,0,(268435456,))]
        expected=[hashlib.sha256(json.dumps(x,separators=(',',':'),ensure_ascii=True).encode()).digest() for x in values]
        with patch.dict(os.environ,MFM_TRACE_HASH_BACKEND='native'):got=source_digests(values)
        self.assertEqual([bytes(x) for x in got],expected)
        for field,value in enumerate(['T',2,3,4,5,1,(48,)]):
            changed=list(values[0]);changed[field]=value
            with patch.dict(os.environ,MFM_TRACE_HASH_BACKEND='native'):other=source_digests([changed])
            self.assertNotEqual(bytes(other[0]),expected[0])

    def test_segment_hash_preserves_coordinate_and_orientation_identity(self):
        reads=np.frombuffer(bytes(range(32))*3,dtype='V32')
        fields=np.array([[0,0,10,20,30,1],[1,5,20,100,115,-1],[2147483647,1,2,3,4,1]],dtype='i8')
        expected=[hashlib.sha256(json.dumps((bytes(r).hex(),*map(int,f)),separators=(',',':')).encode()).digest() for r,f in zip(reads,fields)]
        with patch.dict(os.environ,MFM_TRACE_HASH_BACKEND='native'):actual=segment_digests(reads,fields)
        self.assertEqual([bytes(x) for x in actual],expected)
        for col in range(6):
            changed=fields.copy();changed[0,col]+=1
            with patch.dict(os.environ,MFM_TRACE_HASH_BACKEND='native'):other=segment_digests(reads,changed)
            self.assertNotEqual(bytes(other[0]),expected[0])

    def test_packed_batch_and_invalid_boundaries(self):
        bases=np.array([0x01,0x23,0x40,0x40],dtype='u1');offsets=[0,3,4]
        meta=np.array([[1,30,10,15,0,0,1,5],[0,0,0,1,1,1,0,1]],dtype='i8');ops=[80]
        with patch.dict(os.environ,MFM_TRACE_HASH_BACKEND='native'):a=digests(bases,offsets,meta,ops,packed=True)
        with patch.dict(os.environ,MFM_TRACE_HASH_BACKEND='python'):b=digests(bases,offsets,meta,ops,packed=True)
        np.testing.assert_array_equal(a,b)
        for invalid in [[0,3,9],[0,-1,4]]:
            with self.assertRaises(ValueError),patch.dict(os.environ,MFM_TRACE_HASH_BACKEND='native'):digests(bases,invalid,meta,ops,packed=True)
        bad=bases.copy();bad[0]=0xF0
        with self.assertRaises(ValueError),patch.dict(os.environ,MFM_TRACE_HASH_BACKEND='native'):digests(bad,offsets,meta,ops,packed=True)

    def test_sort_uses_same_68_budget_without_two_thread_cap(self):
        cfg=SimpleNamespace(platform='hifi',threads=23,sort_threads=2,hap_jobs=2,sort_mem='2G')
        env={'MFM_WGS_ALIGN_REQUESTED':'auto','MFM_WGS_SORT_REQUESTED':'auto','MFM_WGS_HAP_REQUESTED':'auto'}
        p=alignment_plan(cfg,{1:3_000_000_000,2:3_000_000_000},68,1000*GIB,env)
        self.assertEqual((p['threads'],p['sort_threads'],p['input_decode_threads']),(17,8,4))
        self.assertEqual(p['reserved_cpus'],68)
        p=alignment_plan(cfg,{1:3_000_000_000},68,1000*GIB,dict(env,MFM_WGS_SORT_REQUESTED='2'))
        self.assertEqual(p['sort_threads'],2)

    def test_gpu_next_hap_fits_while_first_encodes(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{'MFM_WGS_BACKEND':'parabricks'},clear=True):
            phase=SimpleNamespace(path=None,total=68,ceiling=lambda:1000*GIB)
            owner=HaplotypeResources(phase,None,[1,2],local_root=Path(raw))
            env={'MFM_GPU_LOCAL_LEASE_FILE':str(owner.path),'MFM_GPU_HAP_ID':'1','MFM_GPU_STAGE_CPUS':'68'}
            with patch.dict(os.environ,env):
                first=StageBudget();first.output()
            self.assertEqual(first.cpus,16)
            with patch.dict(os.environ,dict(env,MFM_GPU_HAP_ID='2')):
                second=StageBudget();self.assertEqual(second.available_cpus(),52)
                self.assertTrue(second.try_phase('gpu_compute',52,150*GIB))
            row=json.loads(owner.path.read_text());self.assertEqual(row['cpus'],68)
            self.assertEqual(owner.controller(1).reserve('output_qc',max_units=23),16)
            owner.release(2)
            self.assertEqual(owner.controller(1).reserve('output_qc',max_units=23),23)

    def test_existing_idle_loan_is_bounded_by_initial_grant(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();own=root/'leases/own.json'
            _atomic(own,{'cpus':9,'total_cpus':68,'memory_bytes':8*GIB,'cpu_backfill':{'cpus':20},
                'hap_resources_version':2,'hap_resources':{'1':{'cpus':1,'memory_bytes':4*GIB,'phase':'gpu_wait'},'2':{'cpus':8,'memory_bytes':4*GIB,'phase':'output_qc'}}})
            _atomic(root/'leases/other.json',{'cpus':68,'memory_bytes':4*GIB})
            env={'MFM_CPU_POOL_LEASE_FILE':str(own),'MFM_GPU_HAP_ID':'1','MFM_GPU_STAGE_CPUS':'68'}
            with patch.dict(os.environ,env,clear=True):
                budget=StageBudget()
                with patch('integrate_wgs_rna.realign.idle_resources.loan_capacity',return_value=(68,500*GIB,{})):
                    self.assertEqual(budget.available_cpus(),12)
                with patch('integrate_wgs_rna.realign.idle_resources.loan_capacity',return_value=None):
                    self.assertEqual(budget.available_cpus(),0)

    def test_parallel_cram_preserves_all_reads_fields_order_and_tags(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);ref=root/'ref.fa';ref.write_text('>chr1\n'+'ACGT'*500+'\n>chr2\n'+'ACGT'*500+'\n');pysam.faidx(str(ref))
            bam=root/'input.bam'
            header={'HD':{'SO':'coordinate'},'SQ':[{'SN':c,'LN':2000} for c in ['chr1','chr2']],'RG':[{'ID':'old','SM':'old'}]}
            with pysam.AlignmentFile(str(bam),'wb',header=header) as f:
                for rid in range(2):
                    for i in range(20):
                        r=pysam.AlignedSegment();r.query_name=f'read{rid}_{i}';r.query_sequence='ACGT';r.query_qualities=pysam.qualitystring_to_array('IJKL')
                        r.reference_id=rid;r.reference_start=4*(i//2);r.mapping_quality=60;r.cigar=[(0,4)];r.flag=4 if i==19 else 0
                        r.set_tag('HP',rid+1);r.set_tag('PS',42);r.set_tag('OD',i%2);r.set_tag('RG','old');r.set_tag('ZZ','custom');f.write(r)
                r=pysam.AlignedSegment();r.query_name='unplaced';r.query_sequence='NN';r.query_qualities=pysam.qualitystring_to_array('II');r.flag=4;f.write(r)
            outputs=[]
            for cpus in [4,17]:
                target=root/f'{cpus}.cram';details=encode(bam,ref,target,'@RG\\tID:S\\tSM:S','S','preserve_source',cpus=cpus,work=root)
                self.assertEqual(details['backend'],'serial_cram' if cpus==4 else 'native_chromosome_cram')
                if cpus==17:self.assertLessEqual(details['native_cpu_limit'],cpus)
                with pysam.AlignmentFile(str(target),'rc',reference_filename=str(ref)) as f:outputs.append([r.to_string() for r in f])
            self.assertEqual(outputs[0],outputs[1]);self.assertEqual(len(outputs[1]),41)
            self.assertIn('unplaced',outputs[1][-1]);self.assertTrue(all('OD:' not in r for r in outputs[1]))
            for cpus in range(12,69):self.assertLessEqual(plan(cpus,44)['native_cpu_limit'],cpus)

if __name__=='__main__':unittest.main()
