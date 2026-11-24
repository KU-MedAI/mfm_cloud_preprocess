"""QC 계약 회귀: CNV 판정, RNA metadata 호환, config 조기 차단과 경로 분리."""
from __future__ import annotations

import csv
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(os.environ.get("RNA_TEST_ROOT", Path(__file__).resolve().parents[1]))
sys.path[:0] = [str(ROOT), str(ROOT / "tests")]
import rna_profile
from rna_cnv_policy import consume_cnv_row, evaluate_hap_depths
from test_ase_qc import read_qc, run_qc, write_cnv_meta, write_fixture


class QcContractTests(unittest.TestCase):
    def test_cnv_hap_boundary_and_missing(self):
        self.assertEqual(evaluate_hap_depths([0.5, 1.5], 0.5, 1.5), (True, "pass"))
        self.assertEqual(evaluate_hap_depths([0.3, 1.3], 0.5, 1.5), (False, "outside_range"))
        self.assertEqual(evaluate_hap_depths([None, 1.0], 0.5, 1.5), (False, "missing_depth"))
        self.assertEqual(evaluate_hap_depths([None, None], 0.5, 1.5), (False, "missing_depth"))
        for values in ([float("nan"), 1.0], [float("inf"), 1.0], [-1.0, 1.0]):
            with self.assertRaises(ValueError):
                evaluate_hap_depths(values, 0.5, 1.5)

    def test_cnv_producer_flag_and_disabled(self):
        row = {"rd_hap1": "0.3", "rd_hap2": "1.3", "rd": "0.8", "pass_cnv": "0"}
        self.assertEqual(consume_cnv_row(row, enabled=True), (0.8, False, "outside_range"))
        self.assertEqual(consume_cnv_row(None, enabled=True), (None, False, "missing_gene"))
        self.assertEqual(consume_cnv_row(None, enabled=False), (None, True, "disabled"))
        row.update(rd_hap1="NA", pass_cnv="1")  # 구형 TSV의 결측 정상처리를 인정하지 않는다.
        self.assertEqual(consume_cnv_row(row, enabled=True), (0.8, False, "missing_depth"))
        row["pass_cnv"] = "NA"
        with self.assertRaises(ValueError):
            consume_cnv_row(row, enabled=True)

    def test_stage22_cnv_imbalance_missing_and_unaffected_statistics(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_fixture(root, cnv={})
            run_qc(root)  # 필터 비활성 기준선
            before = read_qc(root)
            cnv = root / "gene_rd.tsv"
            with cnv.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            fields = list(rows[0])
            changed = []
            for row in rows:
                if row["gene_id"] == "G1":
                    row.update(rd_hap1="0.3", rd_hap2="1.3", rd="0.8", pass_cnv="0")
                if row["gene_id"] == "G2":
                    continue
                if row["gene_id"] == "G3":
                    row.update(rd_hap1="NA", rd_hap2="NA", rd="NA", pass_cnv="1")
                changed.append(row)
            with cnv.open("w") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(changed)
            write_cnv_meta(root)
            run_qc(root, "--cnv-rd", str(cnv))
            after = read_qc(root)
            self.assertEqual(after["G1"]["pass_cnv"], "0")
            self.assertEqual(after["G1"]["continuous_eligible"], "0")
            self.assertEqual(after["G2"]["cnv_status"], "missing_gene")
            self.assertEqual(after["G3"]["cnv_status"], "missing_depth")
            self.assertEqual(after["G2"]["pass_cnv"], "0")
            for gene in before:
                for key in ("aCount", "bCount", "ratio", "binom_p", "betabinom_p", "fdr_q", "betabinom_fdr_q"):
                    self.assertEqual(before[gene][key], after[gene][key], (gene, key))
            summary = json.loads((root / "ase_qc.json").read_text())
            self.assertEqual(summary["sample_qc"]["genes_cnv_missing_evidence"], 2)
            self.assertFalse(summary["cnv_policy"]["legacy_tsv_without_policy_metadata"])
            self.assertEqual(
                summary["inputs"]["cnv_rd_sha256"],
                json.loads(Path(str(cnv) + ".meta.json").read_text())["output_sha256"],
            )

    def test_stage22_rejects_mismatched_cnv_cutoff(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_fixture(root, cnv={})
            cnv = root / "gene_rd.tsv"
            Path(str(cnv) + ".meta.json").write_text(json.dumps({"min_rd": 0.7}))
            with self.assertRaises(subprocess.CalledProcessError) as caught:
                run_qc(root, "--cnv-rd", str(cnv))
            self.assertIn("cutoff 불일치", caught.exception.stderr)

    def test_coverage_policy_current_and_old_outputs(self):
        current = rna_profile.coverage_metadata()
        rna_profile.validate_coverage_metadata(current)
        self.assertEqual(current["duplicate_policy"], "included")
        self.assertEqual(current["samtools_depth_args"][-2:], ["-g", "1024"])
        old_excluded = {
            "read_filter": "samtools depth: MAPQ255 primary QC-pass; duplicate/supplementary excluded",
            "coverage_policy": "MAPQ 255, primary, QC pass, duplicate 제외, splice N 제외, deletion 포함, mate 겹침 2회",
            "coverage_policy_revision": "2026-10-01_effective_duplicate_excluded",
            "duplicate_policy": "excluded",
            "effective_duplicate_policy": "excluded",
            "legacy_filter_metadata_mismatch": False,
            "samtools_depth_args": ["-a", "-J", "-Q", "255", "-G", "2820"],
        }
        with self.assertRaises(ValueError):
            rna_profile.validate_coverage_metadata(old_excluded)
        for data, key, bad in [(current, "duplicate_policy", "excluded"),
                               (current, "coverage_policy_revision", "unknown"),
                               (current, "samtools_depth_args", ["-a", "-J"]),
                               (current, "legacy_filter_metadata_mismatch", True)]:
            tampered = dict(data, **{key: bad})
            with self.assertRaises(ValueError):
                rna_profile.validate_coverage_metadata(tampered)

    def test_metadata_fix_does_not_invalidate_star_profile_cache(self):
        baseline = os.environ.get("RNA_TEST_BASELINE")
        if not baseline:
            self.skipTest("수정 전 profile 사본이 제공될 때 기준선을 비교합니다")
        spec = importlib.util.spec_from_file_location("baseline_rna_profile", Path(baseline) / "rna_profile.py")
        old = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(old)
        with tempfile.TemporaryDirectory() as temporary:
            gtf = Path(temporary) / "fixture.gtf"
            gtf.write_text("fixture annotation\n")
            self.assertEqual(old.fingerprint(str(gtf), addnl=False), rna_profile.fingerprint(str(gtf), addnl=False))

    def test_stage12_uses_manifest_source_sample_alias(self):
        text = (ROOT / "12_compact_zarr.sh").read_text()
        self.assertIn('SOURCE_SAMPLE_ID=$(sample_source_sample_id "$S")', text)
        self.assertIn('--source-sample-id "$SOURCE_SAMPLE_ID"', text)
        self.assertNotIn('--source-sample-id "$S"', text)

    def config(self, **overrides):
        env = os.environ.copy()
        for key in ("PANEL_BASE", "PANEL_LOCAL_DIR", "PANEL_SOURCE_URL", "WINDOW_DROP_N", "RNA_CNV_FILTER"):
            env.pop(key, None)
        env.update(overrides)
        return subprocess.run(["bash", "-c",
            'source "$1"; printf "%s\\n%s\\n" "$PANEL_LOCAL_DIR" "$PANEL_SOURCE_URL"',
            "test", str(ROOT / "config_v2.sh")], env=env, capture_output=True, text=True)

    def test_config_panel_local_legacy_and_remote(self):
        result = self.config(PANEL_BASE="/tmp/panel")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines()[0], "/tmp/panel")
        result = self.config(PANEL_BASE="https://example.org/panel", PANEL_LOCAL_DIR="/tmp/local")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), ["/tmp/local", "https://example.org/panel"])
        self.assertNotEqual(self.config(PANEL_LOCAL_DIR="https://example.org/panel").returncode, 0)
        self.assertNotEqual(self.config(PANEL_BASE="/tmp/a", PANEL_LOCAL_DIR="/tmp/b").returncode, 0)

    def test_config_accepts_n_preservation_and_rejects_invalid_options(self):
        self.assertEqual(self.config().returncode, 0)
        self.assertEqual(self.config(WINDOW_DROP_N="0").returncode, 0)
        for value in ("2", "bad"):
            result = self.config(WINDOW_DROP_N=value)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("N-window", result.stderr)
        self.assertNotEqual(self.config(RNA_CNV_FILTER="2").returncode, 0)

    def test_local_wrapper_continues_after_failed_sample(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            # /tmp가 noexec인 서버에서도 fixture를 실행할 수 있게 Python을 target으로 쓴다.
            (base / "S1").write_text("raise SystemExit(1)\n")
            (base / "S2").write_text("print('sample-success')\n")
            status = base / "status"
            env = dict(os.environ, RUNNER_TARGET_SCRIPT=sys.executable, RUN_STATUS_DIR=str(status),
                       RUN_LOG_DIR=str(base / "logs"), PYTHON_BIN=sys.executable)
            result = subprocess.run(["bash", "-c", 'printf "S1\\nS2\\n" | xargs -n1 -P1 "$1"',
                "test", str(ROOT / "scripts/run_one_local.sh")], cwd=base, env=env, capture_output=True, text=True)
            logs = "\n".join(path.read_text() for path in sorted((base / "logs").glob("*.log")))
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr + logs)
            self.assertEqual((status / "results/S1.result.tsv").read_text().split("\t")[1], "failed")
            self.assertEqual((status / "results/S2.result.tsv").read_text().split("\t")[1], "success")


if __name__ == "__main__":
    unittest.main()
