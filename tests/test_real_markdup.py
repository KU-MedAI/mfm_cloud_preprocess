"""실제Picard의직렬/동시HAP결과와CPU·heap계획을대조함."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch
import json,os,shutil,subprocess,sys,tempfile,unittest
import pysam
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from integrate_wgs_rna.markdup_resources import markdup_plan,markdup_memory_bytes
STAGE='17' if (ROOT/'09_personal_fasta.sh').exists() else '15'


class MarkdupResourceTest(unittest.TestCase):
    def test_cpu_lease_and_ram_budget(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);bam=root/'reads.bam';bam.write_bytes(b'fixture');lease=root/'lease.json'
            lease.write_text(json.dumps({'cpus':2,'memory_bytes':markdup_memory_bytes(1)}))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_CPU_BUDGET':'72'},clear=True):
                p=markdup_plan(1,bam);self.assertEqual(p['java_active_processors'],2);self.assertEqual(p['htslib_aux_threads'],1)
                self.assertLessEqual(p['estimated_peak_rss_bytes'],p['ram_budget_bytes'])
                with self.assertRaisesRegex(ValueError,'RAM'):markdup_plan(2,bam)

    def test_직접실행도cgroup_CPU상한을지킴(self):
        with tempfile.TemporaryDirectory() as raw:
            bam=Path(raw)/'input.bam';bam.write_bytes(b'input')
            with patch.dict(os.environ,{},clear=True),patch('integrate_wgs_rna.markdup_resources.cgroup_cpu_limit',return_value=2):
                plan=markdup_plan(1,bam)
            self.assertEqual(plan['cpu_grant'],2)
            self.assertEqual(plan['htslib_aux_threads'],1)

    @unittest.skipUnless(shutil.which('picard') and shutil.which('samtools'),'Picard/samtools필요')
    def test_actual_picard_parallel_matches_serial(self):
        # /tmp는noexec일수있으므로JNI임시폴더는코드사본의검사폴더에둠.
        with tempfile.TemporaryDirectory(prefix='.tmp_picard_',dir=ROOT/'tests') as raw:
            root=Path(raw);source=root/'source.bam'
            header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}],
                    'RG':[{'ID':'S','SM':'S','LB':'L','PL':'ILLUMINA','PU':'unit'}]}
            records=[]
            for n,start in enumerate([50,50,300]):
                for flag,pos,mate in [(99,start,start+100),(147,start+100,start)]:
                    r=pysam.AlignedSegment();r.query_name=f'INST:1:FLOW:1:1101:{1000+n}:1000';r.query_sequence='A'*40;r.flag=flag
                    r.reference_id=0;r.reference_start=pos;r.mapping_quality=255;r.cigar=[(0,40)]
                    r.next_reference_id=0;r.next_reference_start=mate;r.template_length=140 if flag==99 else -140
                    r.query_qualities=pysam.qualitystring_to_array('I'*40);r.set_tag('RG','S');records.append(r)
            with pysam.AlignmentFile(str(source),'wb',header=header) as out:
                for r in sorted(records,key=lambda x:x.reference_start):out.write(r)
            env=dict(os.environ,MFM_CPU_BUDGET='4',MFM_SAMPLE_CPU_BUDGET='2',MFM_RAM_LIMIT_BYTES=str(16*1024**3),PYTHONDONTWRITEBYTECODE='1')
            env.pop('MFM_CPU_POOL_LEASE_FILE',None)
            def one(mode,hap):
                d=root/mode/f'HAP{hap}';d.mkdir(parents=True);bam=d/'star.bam';shutil.copyfile(source,bam)
                align=d/'align.json';align.write_text(json.dumps({'align_dir':str(d),'star_bam':str(bam)}))
                cmd=[sys.executable,str(ROOT/f'{STAGE}_mark_duplicates.py'),'--sample','S','--hap',str(hap),'--align-json',str(align),
                     '--out-dir',str(d/'out'),'--work-dir',str(d/'work'),'--log-dir',str(d/'logs'),'--picard',shutil.which('picard'),
                     '--samtools',shutil.which('samtools'),'--heap-gb','1']
                (d/'work').mkdir();result=subprocess.run(cmd,env=env,capture_output=True,text=True,timeout=60)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                metadata=json.loads((d/f'work/hap{hap}.markdup.json').read_text())
                with pysam.AlignmentFile(metadata['marked_bam'],'rb') as handle:
                    rows=[(r.query_name,r.flag,r.reference_start,r.cigarstring,tuple(r.tags)) for r in handle.fetch(until_eof=True)]
                self.assertEqual(metadata['resource_execution']['htslib_aux_threads'],1)
                self.assertEqual(metadata['resource_execution']['java_active_processors'],2)
                return rows,{k:metadata[k] for k in ['primary_duplicate_records','primary_records','read_pair_duplicates','percent_duplication']}
            serial=[one('serial',h) for h in [1,2]]
            with ThreadPoolExecutor(max_workers=2) as pool:parallel=list(pool.map(lambda h:one('parallel',h),[1,2]))
            self.assertEqual(serial,parallel);self.assertEqual(serial[0][1]['primary_duplicate_records'],2)


if __name__=='__main__':unittest.main()
