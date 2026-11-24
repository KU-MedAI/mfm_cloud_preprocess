"""실패 종료 코드·stderr·반복 시도별 로그와 JSONL 연결을 검사한다."""
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class StageFailureLogsTest(unittest.TestCase):
    def test_missing_wgs_cnv_records_stage22_and_next_sample_can_run(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "scripts").mkdir()
            for name in ("record_stage_failure.py", "stage_provenance.py"):
                shutil.copy2(ROOT / "scripts" / name, directory / "scripts" / name)
            exons = directory / "exons.tsv"
            exons.write_text("g1\tchr1\t0\t100\n")
            final = directory / "final"
            script = f'''
set -euo pipefail
source "{ROOT / 'lib/stage_runner.sh'}"
sample_final_root() {{ printf '%s\\n' "{final}"; }}
log() {{ printf '%s\\n' "$*" >&2; }}
pipeline_write_sample_state() {{ :; }}
pipeline_record_stage() {{ :; }}
PIPELINE_SAMPLE_ID=S1
PIPELINE_STAGE_LABEL=rna
PIPELINE_START_AT=22
PIPELINE_STOP_AFTER=
PYTHON="{os.sys.executable}"
RUN_ID=test
rc=0
pipeline_run_command 22 22_cnv_filter.py "$PYTHON" "{ROOT / '22_cnv_filter.py'}" \\
  --exons "{exons}" --bam "{directory / 'missing.cram'}" --out "{directory / 'gene_rd.tsv'}" || rc=$?
[[ "$rc" != 0 ]] || exit 80
PIPELINE_SAMPLE_ID=S2
pipeline_run_command 22 next_sample bash -c 'echo next-sample-success'
'''
            result = subprocess.run(["bash", "-c", script], cwd=directory, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            failures = [json.loads(line) for line in (final / "S1/metadata/stage_failures.jsonl").read_text().splitlines()]
            self.assertEqual(failures[0]["stage"], "22")
            self.assertIn("WGS 정렬 없음", failures[0]["log_tail"])
            events = [json.loads(line) for line in (final / "S2/metadata/execution_provenance.jsonl").read_text().splitlines()]
            self.assertEqual(events[-1]["status"], "success")
            self.assertFalse((directory / "gene_rd.tsv").exists())

    def test_exit_code_logs_and_retries_are_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "scripts").mkdir()
            for name in ("record_stage_failure.py", "stage_provenance.py"):
                shutil.copy2(ROOT / "scripts" / name, directory / "scripts" / name)
            final = directory / "final"
            script = f'''
set -euo pipefail
source "{ROOT / 'lib/stage_runner.sh'}"
sample_final_root() {{ printf '%s\\n' "{final}"; }}
log() {{ printf '%s\\n' "$*" >&2; }}
pipeline_write_sample_state() {{ :; }}
pipeline_record_stage() {{ :; }}
PIPELINE_SAMPLE_ID=S1
PIPELINE_STAGE_LABEL=rna
PIPELINE_START_AT=18
PIPELINE_STOP_AFTER=
PYTHON="{os.sys.executable}"
RUN_ID=test
for attempt in 1 2; do
  rc=0
  pipeline_run_command 18 intentional_failure bash -c 'echo stdout; echo "ValueError: fixture failure" >&2; exit 23' || rc=$?
  [[ "$rc" == 23 ]] || exit 80
done
'''
            result = subprocess.run(["bash", "-c", script], cwd=directory, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            sample = final / "S1"
            failures = [json.loads(line) for line in (sample / "metadata/stage_failures.jsonl").read_text().splitlines()]
            self.assertEqual(len(failures), 2)
            self.assertNotEqual(failures[0]["log_path"], failures[1]["log_path"])
            for failure in failures:
                self.assertEqual(failure["exit_code"], 23)
                self.assertEqual(failure["diagnosis_status"], "unclassified")
                self.assertIn("ValueError: fixture failure", failure["error_signature"])
                self.assertIn("stdout", Path(failure["log_path"]).read_text())
            provenance = [json.loads(line) for line in (sample / "metadata/execution_provenance.jsonl").read_text().splitlines()]
            finish = [event for event in provenance if event["event"] == "finish"]
            self.assertEqual([event["exit_code"] for event in finish], [23, 23])
            self.assertEqual([event["log_path"] for event in finish], [event["log_path"] for event in failures])


if __name__ == "__main__":
    unittest.main()
