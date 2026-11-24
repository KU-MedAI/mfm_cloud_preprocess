"""독립 분기 선행조건·실패·번호 재개와 CPU 상한을 검증함."""
from pathlib import Path
import json,os,subprocess,sys,tempfile,time,unittest
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=r'''set -euo pipefail
source "$ROOT/lib/hap_stage_group.sh"
source "$ROOT/lib/rna_dag.sh"
W=$TEMP;F=$TEMP/final
mkdir -p "$W/rna" "$F/qc"
PIPELINE_START_AT=$START;PIPELINE_STOP_AFTER=
RNA_STAR_INDEX_BENCHMARK=0
RNA_STAR_INDEX_HAP_JOBS=2;RNA_MARKDUP_HAP_JOBS=2;RNA_ALIGN_HAP_JOBS=2;RNA_COVERAGE_HAP_JOBS=2
export MFM_CPU_BUDGET=16
pipeline_write_sample_state() { :; }
pipeline_record_stage() { :; }
pipeline_collect_statistics() { :; }
event() { printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "${3:-0}" "${MFM_CPU_POOL_STAGE_MAX:-0}" "$(date +%s.%N)" >> "$TEMP/events"; }
rna_star_index_hap() { event start index "$1"; if [[ $1 == 2 ]]; then sleep .5; else sleep .05; fi; event end index "$1"; }
rna_star_align_hap() { event start align "$1"; [[ ${TEST_FAIL:-0} != 1 || $1 != 1 ]] || exit 23; sleep .05; event end align "$1"; }
rna_markdup_hap() { event start mark "$1"; sleep .05; event end mark "$1"; }
rna_coverage_hap() { event start coverage "$1"; sleep .5; event end coverage "$1"; }
rna_qc_hap() { event start qc "$1"; event end qc "$1"; }
rna_variants() { event start variants; sleep .1; event end variants; }
rna_cnv() { event start cnv; sleep .15; event end cnv; }
rna_allelic() { event start ase; sleep .15; event end ase; }
rna_dag_enabled
pipeline_rna_dag
event end all
'''


class DagTests(unittest.TestCase):
    def test_fixed_sample_indexes_use_returned_evidence_cpu(self):
        script=SCRIPT.replace('export MFM_CPU_BUDGET=16','export MFM_CPU_BUDGET=16 MFM_CPU_ALLOCATION_MODE=fixed_sample')
        result,folder=self.execute(script=script)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        events=[line.split('\t') for line in (folder/'events').read_text().splitlines()]
        ends={name:float(stamp) for kind,name,h,cap,stamp in events if kind=='end' and name in {'variants','cnv'}}
        for kind,name,h,cap,stamp in events:
            if kind=='start' and name=='index':
                self.assertEqual(int(cap),8)
                self.assertGreaterEqual(float(stamp),max(ends.values()))
        active={}
        for kind,name,h,cap,stamp in sorted(events,key=lambda row:float(row[4])):
            if name=='all':continue
            if kind=='start':active[name,h]=int(cap)
            else:active.pop((name,h),None)
            self.assertLessEqual(sum(active.values()),16)

    def test_fixed_sample_failed_evidence_prevents_index_start(self):
        script=SCRIPT.replace('export MFM_CPU_BUDGET=16','export MFM_CPU_BUDGET=16 MFM_CPU_ALLOCATION_MODE=fixed_sample')
        script=script.replace('sleep .15; event end cnv','exit 9; event end cnv')
        result,folder=self.execute(script=script)
        self.assertNotEqual(result.returncode,0)
        self.assertNotIn('start\tindex',(folder/'events').read_text())

    def execute(self,script=SCRIPT,fail=False):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);folder=Path(self.temp.name)
        run=subprocess.run(['bash','-c',script],env={**os.environ,'ROOT':str(ROOT),'TEMP':str(folder),'PYTHON':sys.executable,'TEST_FAIL':str(int(fail)),'START':'13' if (ROOT/'09_personal_fasta.sh').exists() else '11'},capture_output=True,text=True,timeout=15)
        return run,folder

    def test_dependencies_overlap_and_concurrent_cpu_caps(self):
        result,folder=self.execute();self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        events=[line.split('\t') for line in (folder/'events').read_text().splitlines()]
        def stamp(kind,name,h='0'):return float(next(x[4] for x in events if x[:3]==[kind,name,h]))
        self.assertLess(stamp('start','align','1'),stamp('end','index','2'))
        self.assertLess(stamp('start','variants'),stamp('end','index','2'))
        for h in ['1','2']:self.assertLess(stamp('end','mark',h),stamp('start','ase'))
        self.assertLess(stamp('end','variants'),stamp('start','ase'))
        self.assertLess(stamp('start','ase'),stamp('end','coverage','2'))
        timeline=sorted(events,key=lambda x:float(x[4]));active={}
        for kind,name,h,cap,_ in timeline:
            if name=='all':continue
            key=(name,h)
            if kind=='start':active[key]=int(cap)
            else:active.pop(key)
            self.assertLessEqual(sum(active.values()),16,(active,timeline))
        report=json.loads((folder/'final/qc/rna_dag_execution.json').read_text());self.assertEqual(report['status'],'pass')

    def test_failed_hap_unblocks_ase_and_preserves_valid_sibling(self):
        result,folder=self.execute(fail=True)
        self.assertNotEqual(result.returncode,0)
        events=(folder/'events').read_text();self.assertNotIn('start\tase',events);self.assertIn('end\tqc\t2',events)
        report=json.loads((folder/'final/qc/rna_dag_execution.json').read_text());self.assertEqual(report['status'],'failed')

    def test_signal_terminates_its_descendants(self):
        with tempfile.TemporaryDirectory() as raw:
            folder=Path(raw)
            script=SCRIPT.replace('sleep .5','sleep 30')
            process=subprocess.Popen(['bash','-c',script],env={**os.environ,'ROOT':str(ROOT),'TEMP':str(folder),'PYTHON':sys.executable,'TEST_FAIL':'0','START':'13' if (ROOT/'09_personal_fasta.sh').exists() else '11'},stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True)
            deadline=time.monotonic()+5
            while not (folder/'events').exists() and time.monotonic()<deadline:time.sleep(.02)
            descendants=set()
            for name in os.listdir('/proc'):
                if not name.isdigit():continue
                try:
                    fields=(Path('/proc')/name/'stat').read_text().rsplit(')',1)[1].split()
                    if int(fields[3])==process.pid:descendants.add(int(name))
                except (OSError,ValueError,IndexError):pass
            process.terminate();out,err=process.communicate(timeout=10)
            self.assertEqual(process.returncode,143,out+err)
            for pid in descendants:
                path=Path('/proc')/str(pid)/'stat'
                if path.exists():self.assertIn(path.read_text().rsplit(')',1)[1].split()[0],['Z','X'])

    def test_explicit_stop_restart_and_hap_jobs_disable_overlap(self):
        script=SCRIPT[:SCRIPT.index('rna_dag_enabled\n')]+r'''
for mode in stop restart jobs; do
  PIPELINE_STOP_AFTER=;PIPELINE_START_AT=13;RNA_ALIGN_HAP_JOBS=2
  case "$mode" in stop) PIPELINE_STOP_AFTER=18;; restart) PIPELINE_START_AT=21;; jobs) RNA_ALIGN_HAP_JOBS=1;; esac
  if rna_dag_enabled; then exit 99; fi
done
'''
        result,_=self.execute(script);self.assertEqual(result.returncode,0,result.stdout+result.stderr)


if __name__=='__main__':unittest.main()
