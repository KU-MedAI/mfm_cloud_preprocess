"""실제 index 단계 CLI가 도구 실행 시각·thread를 기록하고 cache를 재사용하는지 검사함."""
from pathlib import Path
import json,os,subprocess,sys,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1]


class IndexExecution(unittest.TestCase):
    def test_index_cli_records_actual_command_span_and_reuses_cache(self):
        # /tmp는 noexec일 수 있어 실행 fixture는 workspace에 둠.
        with tempfile.TemporaryDirectory(prefix='.star_exec_',dir=ROOT.parent) as raw:
            base=Path(raw);hap=base/'hap';hap.mkdir();(hap/'chr1.fa').write_text('>chr1\n'+'A'*200+'\n')
            gtf=base/'genes.gtf';gtf.write_text('chr1\tx\texon\t1\t100\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
            fq=base/'fastq.json';fq.write_text(json.dumps({'max_read_length':100}))
            star=base/'fake_star';star.write_text('#!'+sys.executable+'\n'+'''import sys,time
from pathlib import Path
if '--version' in sys.argv:print('fixture_STAR');raise SystemExit(0)
p=Path(sys.argv[sys.argv.index('--genomeDir')+1]);time.sleep(.02)
for name in ['Genome','SA','SAindex','chrName.txt','chrLength.txt','chrStart.txt','genomeParameters.txt']:(p/name).write_text(name)
''');star.chmod(0o755)
            script=ROOT/('15_star_index.py' if (ROOT/'15_star_index.py').exists() else '13_star_index.py')
            command=[sys.executable,'-B',str(script),'--sample','S','--hap','1','--hap-fasta-dir',str(hap),'--contigs','chr1','--gtf',str(gtf),'--fastq-qc',str(fq),'--work-dir',str(base/'work'),'--log-dir',str(base/'logs'),'--cache-tag','fixture','--sjdb-overhang','99','--threads','4','--star',str(star),'--samtools','samtools']
            first=subprocess.run(command,capture_output=True,text=True);self.assertEqual(first.returncode,0,first.stdout+first.stderr)
            p=base/'work/hap1.index.json';record=json.loads(p.read_text());span=record['index_execution']
            self.assertGreater(span['finished_at'],span['started_at']);self.assertEqual(span['threads'],4)
            second=subprocess.run(command,capture_output=True,text=True);self.assertEqual(second.returncode,0,second.stdout+second.stderr)
            self.assertEqual(json.loads(p.read_text()),record)


if __name__=='__main__':unittest.main()
