"""fragment 분할·전체 count/부가통계·STAR 길이 항의 동일성을 검증함."""
from pathlib import Path
import collections,copy,importlib.util,math,os,random,sys,unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
NUMBER='21' if (ROOT/'21_allelic_counts.py').exists() else '19'
spec=importlib.util.spec_from_file_location('allelic_fragment_contract',ROOT/f'{NUMBER}_allelic_counts.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
from integrate_wgs_rna.allelic_parallel import assign_partitioned


def slow_assign(A,B,ex1,ex2,both,evidence,keys=None):
    import time
    result=mod.assign_fragments(A,B,ex1,ex2,both,evidence,keys=keys)
    if keys:
        end=time.process_time()+.04
        while time.process_time()<end:pass
    return result


class Evidence:
    def validated(self,h,slot):return True
    def unvalidated(self,h,slot):return False


class FragmentTests(unittest.TestCase):
    def test_complete_counts_stats_and_split_buckets_match_serial(self):
        rng=random.Random(42);A={};B={}
        for i in range(8000):
            key=f'fragment{i}';a=rng.randrange(10);b=rng.randrange(10)
            if i%11:A[key]=[(a,'chr1',[(10,40)]),None]
            if i%13:B[key]=[(b,'chr1',[(10,40)]),None]
        ex=mod.build_exon_index({'chr1':[(0,100,'G1'),(5,80,'G2')]})
        expected=mod.assign_fragments(A,B,ex,ex,False,Evidence())
        actual,plan=assign_partitioned(mod.assign_fragments,A,B,ex,ex,False,Evidence(),workers=3,batch_size=500,benchmark=True)
        self.assertEqual(actual,expected);self.assertTrue(plan['serial_equal'])
        self.assertEqual(plan['mode'],'fork_fragments')
        self.assertGreater(sum(actual[0]['a'].values())+sum(actual[0]['b'].values()),0)
        for only in [False,True]:
            expected=mod.assign_fragments(A,B,ex,ex,only,Evidence())
            actual,report=assign_partitioned(mod.assign_fragments,A,B,ex,ex,only,Evidence(),workers=1)
            self.assertEqual(actual,expected)

    def test_length_term_matches_original_for_edge_spans_and_empty_blocks(self):
        # 수정 전 구현을 그대로 재현하여 부동소수 경계·단일 mate·다른 contig를 대조함.
        def legacy(table,scale):
            if scale==0:return 0
            n=0
            for slot in table.values():
                mates=[m for m in slot if m is not None]
                if not mates or len({m[1] for m in mates})!=1:continue
                lo=min(b[0][0] for m in mates for b in [m[2]] if b)
                hi=max(b[-1][1] for m in mates for b in [m[2]] if b)
                term=mod.star_length_term(hi-lo,scale)
                if term:
                    for i,m in enumerate(slot):
                        if m is not None:slot[i]=(m[0]-term,m[1],m[2])
                    n+=1
            return n
        cases={}
        for i,span in enumerate([1,3,4,5,63,64,65,1023,1024,1025,16384,262144,262145]):
            cases[str(i)]=[(20,'chr1',[(0,span)]),(15,'chr1',[(span//2,span)])]
        cases.update(single=[None,(12,'chr1',[(5,15)])],different=[(5,'chr1',[(0,5)]),(5,'chr2',[(0,10)])],empty=[None,None],one_empty=[(1,'chr1',[]),(2,'chr1',[(4,8)])])
        for scale in [-.25,0,.25,-1]:
            a,b=copy.deepcopy(cases),copy.deepcopy(cases)
            self.assertEqual(mod.remove_length_term(a,scale),legacy(b,scale));self.assertEqual(a,b)
        for fn in [legacy,mod.remove_length_term]:
            with self.assertRaises(ValueError):fn({'bad':[(1,'chr1',[]),None]},-.25)


    def test_shared_pool_releases_cpu_and_worker_ram_after_join(self):
        import json,tempfile
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as raw:
            pool=Path(raw);(pool/'leases').mkdir();(pool/'waiters').mkdir();lease=pool/'leases/own.json'
            lease.write_text(json.dumps({'cpus':4,'total_cpus':4,'memory_bytes':64*1024**3}))
            A={str(i):[(2,'chr1',[(0,10)]),None] for i in range(3000)}
            B={key:[(1,'chr1',[(0,10)]),None] for key in A}
            ex=mod.build_exon_index({'chr1':[(0,50,'G')]})
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_RAM_LIMIT_BYTES':'1000000000000'}):
                result,plan=assign_partitioned(mod.assign_fragments,A,B,ex,ex,False,Evidence(),workers=3,batch_size=500)
            row=json.loads(lease.read_text());self.assertEqual(row['cpus'],1)
            self.assertEqual(row['memory_bytes'],plan['parent_reserved_bytes'])
            self.assertEqual(sum(result[0]['a'].values()),3000)

    def test_execution_overlap_summary_distinguishes_serial_and_parallel(self):
        from integrate_wgs_rna.allelic_parallel import execution_summary
        a={'pid':1,'started':0.,'finished':5.,'cpu_seconds':4.}
        b={'pid':2,'started':2.,'finished':6.,'cpu_seconds':3.}
        result=execution_summary([a,b]);self.assertEqual(result['peak_active_workers'],2);self.assertTrue(result['parallel_observed']);self.assertEqual(result['worker_cpu_seconds'],7.)
        b.update(started=5.,finished=9.)
        self.assertFalse(execution_summary([a,b])['parallel_observed'])

    def test_real_workers_report_cpu_and_actual_overlap(self):
        A={str(i):[(2,'chr1',[(0,10)]),None] for i in range(4000)}
        B={key:[(1,'chr1',[(0,10)]),None] for key in A}
        ex=mod.build_exon_index({'chr1':[(0,50,'G')]})
        value,plan=assign_partitioned(slow_assign,A,B,ex,ex,False,Evidence(),workers=2,batch_size=500,benchmark=True)
        self.assertTrue(plan['serial_equal']);self.assertTrue(plan['parallel_observed']);self.assertGreater(plan['children_cpu_seconds'],0)


if __name__=='__main__':unittest.main()
