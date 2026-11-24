"""실물 index 비교가 두 HAP 완료 뒤 STOP을 평가하고 캐시를 보존하는지 검증함."""
from pathlib import Path
import json,os,subprocess,sys,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1]


class StageTests(unittest.TestCase):
    def test_serial_parallel_benchmark_and_stop_barrier(self):
        with tempfile.TemporaryDirectory() as raw:
            temp=Path(raw);script=temp/'run.sh'
            (temp/'writer.py').write_text('''import json,sys
from pathlib import Path
h,mode,root=sys.argv[1:];p=Path(root);idx=p/(mode+h);idx.mkdir(exist_ok=True)
for name in ['Genome','SA','SAindex','chrName.txt','chrLength.txt','chrStart.txt']:(idx/name).write_text(name)
(p/'rna'/('hap'+h+'.index.json')).write_text(json.dumps({'index_dir':str(idx)}))
''')
            script.write_text('''set -euo pipefail
source "$ROOT/lib/rna_index_stage.sh"
W=$TEMP
F=$TEMP/final
PIPELINE_DIR=$ROOT
PIPELINE_START_AT=13
PIPELINE_STOP_AFTER=15
RNA_STAR_INDEX_BENCHMARK=1
RNA_STAR_INDEX_HAP_JOBS=1
MFM_SAMPLE_CPU_BUDGET=8
RNA_INDEX_FORCE_ARGS=()
mkdir -p "$W/rna" "$F/qc"
pipeline_run_hap_group() { printf '%s\n' "$MFM_CPU_POOL_STAGE_MAX" >> "$W/caps"; TEST_INDEX_MODE=serial; [[ $3 == 1 ]] || TEST_INDEX_MODE=parallel; printf '%s\\n' "$3" >> "$W/groups"; "$2" 1; "$2" 2; }
pipeline_finish_if_stop() { [[ "$PIPELINE_STOP_AFTER" == 15 ]]; touch "$W/stopped"; }
callback() { mode=$TEST_INDEX_MODE; [[ "${RNA_INDEX_FORCE_ARGS[0]}" == --force ]]; "$PYTHON" "$TEMP/writer.py" "$1" "$mode" "$TEMP"; }
pipeline_run_rna_index_stage 15 callback
[[ ${#RNA_INDEX_FORCE_ARGS[@]} == 0 ]]
''')
            result=subprocess.run(['bash',str(script)],env=dict(os.environ,ROOT=str(ROOT),TEMP=str(temp),PYTHON=sys.executable),capture_output=True,text=True,timeout=20)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertEqual((temp/'groups').read_text().splitlines(),['1','2'])
            self.assertEqual((temp/'caps').read_text().splitlines(),['8','4'])
            self.assertTrue((temp/'stopped').is_file())
            self.assertEqual(json.loads((temp/'final/qc/rna_index_benchmark.json').read_text())['status'],'pass')

    def test_cpu_cap_environment_is_restored_and_failures_propagate(self):
        script='set -euo pipefail\nsource "$ROOT/lib/rna_index_stage.sh"\nMFM_SAMPLE_CPU_BUDGET=8\nMFM_CPU_POOL_STAGE_MAX=6\npipeline_run_hap_group() { [[ "$MFM_CPU_POOL_STAGE_MAX" == 3 ]]; }\npipeline_rna_index_group 15 true 2\n[[ "$MFM_CPU_POOL_STAGE_MAX" == 6 ]]\nunset MFM_CPU_POOL_STAGE_MAX\npipeline_run_hap_group() { return 17; }\nrc=0; pipeline_rna_index_group 15 true 2 || rc=$?\n[[ "$rc" == 17 && -z ${MFM_CPU_POOL_STAGE_MAX+x} ]]\n'
        result=subprocess.run(['bash','-c',script],env=dict(os.environ,ROOT=str(ROOT)),capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)


if __name__=='__main__':unittest.main()
