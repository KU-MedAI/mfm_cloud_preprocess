#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

WORK = Path(__file__).resolve().parent
CODE = WORK if (WORK / "23_ase_qc.py").is_file() else WORK.parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


# 현재 runner가 실행하는 canonical Stage 21/23을 검증한다.
qc = load_module("ase_qc", CODE / "23_ase_qc.py")
counts_module = load_module("allelic_counts", CODE / "21_allelic_counts.py")
cnv_policy = load_module("rna_cnv_policy", CODE / "rna_cnv_policy.py")

META_OK = {
    "schema": "mfm_allelic_counts_meta_v5",
    "script_revision": "2026-10-02_validated_hap_variant_evidence",
    "counting_unit": "fragment",
    "include_duplicates": False,
    "fragments_total": 1000, "fragments_assigned": 100, "fragments_tie": 800,
    "fragments_no_common_mate": 10, "fragments_only_hap1": 5, "fragments_only_hap2": 5,
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_cnv_meta(root: Path, *, min_rd: float = 0.5, max_rd: float = 1.5) -> None:
    path = root / "gene_rd.tsv"
    Path(str(path) + ".meta.json").write_text(json.dumps({
        "schema": "mfm_cnv_rd_policy_v1",
        "policy_revision": cnv_policy.CNV_POLICY_REVISION,
        "min_rd": min_rd,
        "max_rd": max_rd,
        "hap_count": 2,
        "missing_policy": "exclude_from_training",
        "output": str(path),
        "output_sha256": digest(path),
    }))


def write_fixture(root: Path, *, n_filler: int = 400, meta: dict | None = None,
                  cnv: dict[str, float] | None = None,
                  mmap_survives: dict[str, int] | None = None) -> None:
    """G1/G2/G3 + overdispersion 적합에 필요한 균형 발현 filler 유전자."""
    rng = random.Random(11)
    rows = [("G1", 18, 2), ("G2", 10, 10), ("G3", 18, 2)]
    split = {
        "G1": (9, 1, 9, 1),
        "G2": (5, 5, 5, 5),
        "G3": (10, 0, 8, 2),
    }
    for i in range(n_filler):
        total = rng.randrange(20, 80, 2)
        a = sum(1 for _ in range(total) if rng.random() < 0.5)
        gid = f"F{i:04d}"
        rows.append((gid, a, total - a))
        half = total // 2
        a0 = sum(1 for _ in range(half) if rng.random() < 0.5)
        a1 = sum(1 for _ in range(half) if rng.random() < 0.5)
        split[gid] = (a0, half - a0, a1, half - a1)

    with (root / "allelic_counts.tsv").open("w") as fh:
        fh.write("gene_id\taCount\tbCount\ttotalCount\tratio\n")
        for gid, a, b in rows:
            t = a + b
            fh.write(f"{gid}\t{a}\t{b}\t{t}\t{a / t:.6f}\n")

    with (root / "allelic_counts.split_half.tsv").open("w") as fh:
        fh.write("gene_id\tsplit0_aCount\tsplit0_bCount\tsplit0_totalCount\tsplit0_ratio\t"
                 "split1_aCount\tsplit1_bCount\tsplit1_totalCount\tsplit1_ratio\n")
        for gid, _, _ in rows:
            a0, b0, a1, b1 = split[gid]
            t0, t1 = a0 + b0, a1 + b1
            fh.write(f"{gid}\t{a0}\t{b0}\t{t0}\t{a0 / t0 if t0 else 0.5:.6f}\t"
                     f"{a1}\t{b1}\t{t1}\t{a1 / t1 if t1 else 0.5:.6f}\n")

    with (root / "allelic_counts.both_haps.tsv").open("w") as fh:
        fh.write("gene_id\taCount\tbCount\ttotalCount\tratio\n")
        for gid, a, b in rows:
            t = a + b
            fh.write(f"{gid}\t{a}\t{b}\t{t}\t{a / t:.6f}\n")

    with (root / "allelic_counts.mmap_adjusted.tsv").open("w") as fh:
        fh.write("gene_id\taCount\tbCount\taMmap\tbMmap\t"
                 "adj_aCount\tadj_bCount\tadj_ratio\tmmap_survives\n")
        for gid, a, b in rows:
            surv = (mmap_survives or {}).get(gid, 1)
            fh.write(f"{gid}\t{a}\t{b}\t0\t0\t{a}\t{b}\t{a / (a + b):.6f}\t{surv}\n")

    if cnv is not None:
        with (root / "gene_rd.tsv").open("w") as fh:
            fh.write("gene_id\texon_bp\trd_hap1\trd_hap2\trd\tpass_cnv\n")
            for gid, a, b in rows:
                rd = cnv.get(gid, 1.0)
                fh.write(f"{gid}\t1000\t{rd:.4f}\t{rd:.4f}\t{rd:.4f}\t"
                         f"{1 if 0.5 <= rd <= 1.5 else 0}\n")
        write_cnv_meta(root)

    hap_variants = root / "hap_variants.tsv"
    hap_variants.write_text(
        "contig\tstatus\tkind\th1_start\th1_end\th2_start\th2_end\tgenes\n"
        "chr1\tvalidated\tSNV\t100\t101\t100\t101\tG1\n"
        "chr1\tunsupported\tINDEL\t200\t201\t200\t202\tG3\n"
    )
    meta_data = dict(meta or META_OK)
    meta_data.setdefault("hap_variants_sha256", hashlib.sha256(hap_variants.read_bytes()).hexdigest())
    meta_data["output_sha256"] = {
        "counts": digest(root / "allelic_counts.tsv"),
        "split_half": digest(root / "allelic_counts.split_half.tsv"),
        "both_haps": digest(root / "allelic_counts.both_haps.tsv"),
        "mmap_adjusted": digest(root / "allelic_counts.mmap_adjusted.tsv"),
    }
    (root / "allelic_counts.meta.json").write_text(json.dumps(meta_data))


def run_qc(root: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run([
        sys.executable, str(CODE / "23_ase_qc.py"),
        "--counts", str(root / "allelic_counts.tsv"),
        "--hap-variants", str(root / "hap_variants.tsv"),
        "--mmap-adjusted", str(root / "allelic_counts.mmap_adjusted.tsv"),
        "--max-split-ratio-delta", "0.1",
        "--min-eligible-genes", "2",
        "--min-split-correlation", "-1",
        *extra,
    ], check=True, capture_output=True, text=True)


def read_qc(root: Path) -> dict[str, dict[str, str]]:
    with (root / "ase_qc.tsv").open() as fh:
        return {r["gene_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


class TestAseQc(unittest.TestCase):
    def test_exact_binomial_and_bh(self):
        self.assertAlmostEqual(qc.binom_two_sided(10, 0), 2 / 1024)
        self.assertEqual(qc.binom_two_sided(5, 5), 1.0)
        adjusted = qc.bh_fdr([0.01, 0.04, 0.03])
        self.assertEqual([round(x, 3) for x in adjusted], [0.03, 0.04, 0.04])

    def test_split_is_deterministic(self):
        self.assertEqual(counts_module.split_bucket("fragment-1"),
                         counts_module.split_bucket("fragment-1"))
        self.assertIn(counts_module.split_bucket("fragment-2"), (0, 1))

    # ── beta-binomial (AlleleSeq2 alleledb_*) ────────────────────────────────
    def test_betabinomial_reduces_to_binomial_at_rho_zero(self):
        for a, b in [(10, 2), (50, 50), (3, 20), (100, 70)]:
            self.assertAlmostEqual(qc.betabinom_two_sided(a, b, 0.0),
                                   qc.binom_two_sided(a, b), places=12)

    def test_betabinomial_pmf_sums_to_one(self):
        for n, rho in [(10, 0.05), (25, 0.2), (7, 0.5)]:
            self.assertAlmostEqual(sum(qc.betabinom_pmf_series(n, 0.5, rho, n)), 1.0, places=12)

    def test_overdispersion_makes_pvalues_conservative(self):
        p0 = qc.betabinom_two_sided(80, 40, 0.0)
        for rho in (0.05, 0.2, 0.4):
            self.assertGreater(qc.betabinom_two_sided(80, 40, rho), p0)

    def test_overdispersion_fit_recovers_rho(self):
        rng = random.Random(3)
        rows = []
        true_rho = 0.10
        s = (1 - true_rho) / true_rho
        for _ in range(3000):
            n = rng.randrange(20, 200)
            p = rng.betavariate(0.5 * s, 0.5 * s)
            a = sum(1 for _ in range(n) if rng.random() < p)
            rows.append((a, n))
        rho, diag = qc.fit_overdispersion(rows)
        self.assertEqual(diag["status"], "ok")
        self.assertAlmostEqual(rho, true_rho, delta=0.03)

    def test_fit_refuses_small_gene_sets(self):
        rho, diag = qc.fit_overdispersion([(5, 10), (6, 12), (7, 14)])
        self.assertEqual(rho, 0.0)
        self.assertEqual(diag["status"], "insufficient_data")

    # ── label 규칙 ───────────────────────────────────────────────────────────
    def test_continuous_label_does_not_require_fdr(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            write_fixture(root)
            run_qc(root)
            rows = read_qc(root)
            # 연속형은 깊이·안정성만 본다 — 유의성과 무관하다.
            self.assertEqual(rows["G1"]["continuous_eligible"], "1")
            self.assertEqual(rows["G2"]["continuous_eligible"], "1")
            # G3는 split 비율이 흔들려 제외된다.
            self.assertEqual(rows["G3"]["continuous_eligible"], "0")
            # G2는 완전 균형이라 어느 검정으로도 유의하지 않다.
            self.assertEqual(rows["G2"]["binary_ase_eligible"], "0")
            summary = json.loads((root / "ase_qc.json").read_text())
            self.assertTrue(summary["sample_qc"]["sample_training_ready"])
            self.assertEqual(summary["schema"], "mfm_ase_qc_v3")
            self.assertEqual(summary["statistics"]["overdispersion_fit"]["status"], "ok")

    def test_betabinomial_is_the_reported_primary_test(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            write_fixture(root)
            run_qc(root)
            rows = read_qc(root)
            self.assertIn("betabinom_p", rows["G1"])
            self.assertIn("betabinom_fdr_q", rows["G1"])
            self.assertIn("ase_significant_binom", rows["G1"])
            summary = json.loads((root / "ase_qc.json").read_text())
            stats = summary["statistics"]
            self.assertIn("beta-binomial", stats["primary_test"])
            self.assertEqual(stats["overdispersion_fit"]["status"], "ok")
            # overdispersion 보정은 binomial보다 많이 잡을 수 없다.
            self.assertLessEqual(stats["genes_significant_betabinom"],
                                 stats["genes_significant_binom"])

    def test_cnv_region_excludes_gene(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            write_fixture(root, cnv={"G1": 2.4})     # AlleleSeq2 구간 0.5–1.5 밖
            run_qc(root, "--cnv-rd", str(root / "gene_rd.tsv"))
            rows = read_qc(root)
            self.assertEqual(rows["G1"]["pass_cnv"], "0")
            self.assertEqual(rows["G1"]["continuous_eligible"], "0")
            self.assertIn("cnv_region", rows["G1"]["reason_codes"])
            self.assertEqual(rows["G2"]["pass_cnv"], "1")

    def test_mmap_explained_imbalance_is_flagged(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            write_fixture(root, mmap_survives={"G1": 0})
            run_qc(root)
            rows = read_qc(root)
            self.assertEqual(rows["G1"]["pass_mmap_adjustment"], "0")
            self.assertEqual(rows["G1"]["binary_ase_eligible"], "0")
            self.assertIn("mmap_explains_imbalance", rows["G1"]["reason_codes"])

    # ── provenance ───────────────────────────────────────────────────────────
    def test_stale_revision_blocks_sample(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            stale = dict(META_OK, schema="mfm_allelic_counts_meta_v2",
                         script_revision="2026-09-16_fragment_split_half")
            write_fixture(root, meta=stale)
            run_qc(root)
            summary = json.loads((root / "ase_qc.json").read_text())
            self.assertFalse(summary["sample_qc"]["sample_training_ready"])
            self.assertIn("invalid_provenance", summary["sample_qc"]["reason_codes"])

    def test_hap_variant_hash_mismatch_blocks_sample(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            write_fixture(root)
            path = root / "hap_variants.tsv"
            path.write_text(path.read_text().replace("\tG1\n", "\tG2\n"))
            run_qc(root)
            summary = json.loads((root / "ase_qc.json").read_text())
            self.assertFalse(summary["sample_qc"]["sample_training_ready"])
            self.assertIn("invalid_provenance", summary["sample_qc"]["reason_codes"])

    def test_duplicates_included_blocks_sample(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            write_fixture(root, meta=dict(META_OK, include_duplicates=True))
            run_qc(root)
            summary = json.loads((root / "ase_qc.json").read_text())
            self.assertFalse(summary["sample_qc"]["sample_training_ready"])
            self.assertIn("duplicates_included", summary["sample_qc"]["reason_codes"])


if __name__ == "__main__":
    unittest.main()
