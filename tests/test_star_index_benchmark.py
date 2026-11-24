"""STAR 핵심 index 비교·보수 fallback·실제 소형 index 동등성을 검증함."""
from pathlib import Path
import concurrent.futures,json,random,shutil,subprocess,sys,tempfile,time,unittest
ROOT=Path(__file__).resolve().parents[1];HELPER=ROOT/'src/integrate_wgs_rna/rna_index_benchmark.py'


def run_helper(mode,work,seconds):
    return subprocess.run([sys.executable,str(HELPER),mode,'--work-dir',str(work),'--report',str(work/'report.json'),'--seconds',str(seconds)],capture_output=True,text=True,check=True)


class IndexTests(unittest.TestCase):
    def test_core_match_and_mismatch_restores_serial_records(self):
        with tempfile.TemporaryDirectory() as raw:
            work=Path(raw)
            for h in [1,2]:
                directory=work/f'idx{h}';directory.mkdir()
                for name in ['Genome','SA','SAindex','chrName.txt','chrLength.txt','chrStart.txt']:(directory/name).write_text(name)
                (work/f'hap{h}.index.json').write_text(json.dumps({'index_dir':str(directory),'test':'serial','index_execution':{'started_at':(h-1)*10,'finished_at':h*10,'threads':4}}))
            run_helper('snapshot',work,10)
            for h in [1,2]:
                record=json.loads((work/f'hap{h}.index.json').read_text());record['index_execution']={'started_at':20+h,'finished_at':26+h,'threads':2}
                (work/f'hap{h}.index.json').write_text(json.dumps(record))
            run_helper('compare',work,6)
            report=json.loads((work/'report.json').read_text());self.assertEqual(report['recommended_hap_jobs'],2)
            (work/'idx2/Genome').write_text('오염')
            (work/'hap2.index.json').write_text(json.dumps({'index_dir':str(work/'idx2'),'test':'parallel'}))
            run_helper('compare',work,5)
            report=json.loads((work/'report.json').read_text());self.assertFalse(report['core_sha256_equal'])
            self.assertEqual(report['recommended_hap_jobs'],1)
            self.assertEqual(json.loads((work/'hap2.index.json').read_text())['test'],'serial')

    def test_sequential_or_larger_cpu_budget_does_not_select_parallel(self):
        with tempfile.TemporaryDirectory() as raw:
            work=Path(raw)
            def records(starts,threads):
                for h,start in enumerate(starts,1):
                    directory=work/f'idx{h}';directory.mkdir(exist_ok=True)
                    for n in ['Genome','SA','SAindex','chrName.txt','chrLength.txt','chrStart.txt']:(directory/n).write_text(n)
                    (work/f'hap{h}.index.json').write_text(json.dumps({'index_dir':str(directory),'index_execution':{'started_at':start,'finished_at':start+5,'threads':threads}}))
            records([0,10],4);run_helper('snapshot',work,20)
            records([20,25],2);run_helper('compare',work,1)
            report=json.loads((work/'report.json').read_text());self.assertEqual(report['recommended_hap_jobs'],1);self.assertEqual(report['actual_overlap_seconds'],0)
            records([30,31],4);run_helper('compare',work,1)
            report=json.loads((work/'report.json').read_text());self.assertEqual(report['recommended_hap_jobs'],1);self.assertFalse(report['cpu_budget_comparable'])
            records([40,41],2);run_helper('compare',work,1)
            report=json.loads((work/'report.json').read_text());self.assertEqual(report['recommended_hap_jobs'],2);self.assertTrue(report['cpu_budget_comparable'])

    @unittest.skipUnless(shutil.which('STAR'),'STAR 필요')
    def test_real_star_one_and_two_thread_indexes_are_identical(self):
        # 입력100kb·SAindex4, 병렬2CPU·예상RAM2GB 미만의 소형 실제 도구 검사임.
        with tempfile.TemporaryDirectory() as raw:
            work=Path(raw);rng=random.Random(42);fasta=work/'genome.fa';fasta.write_text('>chr1\n'+''.join(rng.choices('ACGT',k=100000))+'\n')
            gtf=work/'genes.gtf';gtf.write_text('chr1\tx\texon\t100\t400\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
            def build(h,mode,threads):
                directory=work/f'{mode}{h}';directory.mkdir()
                command=['STAR','--runMode','genomeGenerate','--runThreadN',str(threads),'--genomeDir',str(directory),'--genomeFastaFiles',str(fasta),'--sjdbGTFfile',str(gtf),'--sjdbOverhang','19','--genomeSAindexNbases','4','--genomeChrBinNbits','10','--outFileNamePrefix',str(directory)+'/']
                started=time.time()
                subprocess.run(command,capture_output=True,check=True)
                (work/f'hap{h}.index.json').write_text(json.dumps({'index_dir':str(directory),'index_execution':{'started_at':started,'finished_at':time.time(),'threads':threads}}))
            for h in [1,2]:build(h,'serial',2)
            run_helper('snapshot',work,2)
            with concurrent.futures.ThreadPoolExecutor(2) as pool:list(pool.map(lambda h:build(h,'parallel',1),[1,2]))
            run_helper('compare',work,1)
            self.assertTrue(json.loads((work/'report.json').read_text())['core_sha256_equal'])


if __name__=='__main__':unittest.main()
