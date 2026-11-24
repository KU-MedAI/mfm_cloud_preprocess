"""실제공유풀을통해HAP병렬·barrier·stop·실패·신호정리를검증함."""
from pathlib import Path
import csv,json,os,signal,subprocess,sys,tempfile,time,unittest
ROOT=Path(__file__).resolve().parents[1]
STAGE='17' if (ROOT/'09_personal_fasta.sh').exists() else '15'


class GroupTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='mfm_group_');self.base=Path(self.temp.name)
        (self.base/'task.py').write_text('''import json,os,sys,time
from pathlib import Path
h,root,fail,duration=sys.argv[1:];p=Path(root)/('hap'+h+'.json')
d={'start':time.time(),'pid':os.getpid(),'grant':int(os.environ.get('MFM_SAMPLE_CPU_BUDGET','0'))}
p.write_text(json.dumps(d));time.sleep(float(duration));d['end']=time.time();p.write_text(json.dumps(d))
raise SystemExit(7 if h==fail else 0)
''')
        self.script=self.base/'run.sh'
        self.script.write_text('''set -euo pipefail
cd "$ROOT"
source lib/state.sh
source lib/stage_runner.sh
source lib/hap_stage_group.sh
log() { printf '%s\\n' "$*"; }
warn() { log "$*"; }
sample_final_root() { printf '%s/final\\n' "$TEST_ROOT"; }
PIPELINE_DIR=$ROOT
PIPELINE_SAMPLE_ID=S
PIPELINE_STAGE_LABEL=RNA
PIPELINE_START_AT=${TEST_START:-$STAGE}
PIPELINE_STOP_AFTER=${TEST_STOP:-}
CURRENT_STAGE=init
LAST_SUCCESS_STAGE=previous
STATE_FILE=$TEST_ROOT/state.tsv
STAGE_FILE=$TEST_ROOT/stages.tsv
PIPELINE_TOOL_VERSIONS_RECORDED=1
pipeline_install_sample_trap
callback() {
  pipeline_run_command "$STAGE" "markdup:HAP$1" "$PYTHON" "$TEST_ROOT/task.py" "$1" "$TEST_ROOT" "${TEST_FAIL:-0}" "${TEST_DURATION:-0.2}"
}
pipeline_run_hap_group "$STAGE" callback "$TEST_JOBS"
printf '%s\\n' next > "$TEST_ROOT/next"
pipeline_complete_sample
''')

    def tearDown(self):self.temp.cleanup()

    def env(self,**kwargs):
        env=dict(os.environ,ROOT=str(ROOT),TEST_ROOT=str(self.base),STAGE=STAGE,TEST_JOBS='2',PYTHON=sys.executable,
                 MFM_CPU_BUDGET='8',MFM_CPU_POOL_DIR=str(self.base/'pool'),MFM_RAM_LIMIT_BYTES=str(64*1024**3),
                 RNA_MARKDUP_HEAP_GB='1',MFM_CPU_POOL_LAYOUT='original' if STAGE=='17' else 'legacy')
        env.pop('MFM_CPU_POOL_LEASE_FILE',None);env.update(kwargs);return env

    def run_case(self,**kwargs):
        out=subprocess.run(['bash',str(self.script)],env=self.env(**kwargs),capture_output=True,text=True,timeout=30)
        return out

    def records(self):return [json.loads((self.base/f'hap{h}.json').read_text()) for h in [1,2]]

    def test_parallel_barrier_and_serial_option(self):
        result=self.run_case();self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        a,b=self.records();self.assertLess(max(a['start'],b['start']),min(a['end'],b['end']))
        self.assertEqual((a['grant'],b['grant']),(4,4))
        state=list(csv.DictReader((self.base/'state.tsv').open(),delimiter='\t'))[0]
        self.assertEqual(state['status'],'success');self.assertEqual(state['last_success_stage'],STAGE)
        self.assertFalse(list((self.base/'pool/leases').glob('*.json')))
        result=self.run_case(TEST_JOBS='1');self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        a,b=self.records();self.assertLessEqual(a['end'],b['start'])

    def test_stop_waits_for_both_haps(self):
        result=self.run_case(TEST_STOP=STAGE);self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertTrue(all('end' in r for r in self.records()));self.assertFalse((self.base/'next').exists())
        rows=list(csv.DictReader((self.base/'final/S/metadata/stage_runtime.tsv').open(),delimiter='\t'))
        self.assertEqual(len(rows),2);self.assertTrue(all(x['status']=='success' for x in rows))

    def test_failure_keeps_sibling_success_and_blocks_next(self):
        result=self.run_case(TEST_FAIL='2');self.assertEqual(result.returncode,7,result.stdout+result.stderr)
        self.assertFalse((self.base/'next').exists());self.assertTrue(all('end' in r for r in self.records()))
        state=list(csv.DictReader((self.base/'state.tsv').open(),delimiter='\t'))[0]
        self.assertEqual(state['status'],'failed');self.assertEqual(state['last_success_stage'],'previous')
        rows=list(csv.DictReader((self.base/'final/S/metadata/stage_runtime.tsv').open(),delimiter='\t'))
        self.assertEqual(sorted(x['status'] for x in rows),['failed','success'])
        self.assertFalse(list((self.base/'pool/leases').glob('*.json')))

    def test_start_skip_and_signal_child_cleanup(self):
        result=self.run_case(TEST_START=str(int(STAGE)+1));self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertFalse((self.base/'hap1.json').exists())
        (self.base/'next').unlink()
        with (self.base/'signal.log').open('w') as log:
            process=subprocess.Popen(['bash',str(self.script)],env=self.env(TEST_DURATION='20'),stdout=log,stderr=subprocess.STDOUT)
            try:
                deadline=time.monotonic()+15
                while not all((self.base/f'hap{h}.json').exists() for h in [1,2]):
                    if process.poll() is not None or time.monotonic()>deadline:raise AssertionError('자식 시작 실패')
                    time.sleep(.05)
                process.send_signal(signal.SIGTERM);self.assertEqual(process.wait(timeout=15),143)
                self.assertFalse((self.base/'next').exists())
                for record in self.records():
                    with self.assertRaises(ProcessLookupError):os.kill(record['pid'],0)
                self.assertFalse(list((self.base/'pool/leases').glob('*.json')))
            finally:
                if process.poll() is None:process.kill();process.wait()


if __name__=='__main__':unittest.main()
