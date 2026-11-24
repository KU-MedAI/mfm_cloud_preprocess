"""N·blacklist 보존부터 4bit 복원·RNA 좌표·최종 mapping까지 검증한다."""
import csv
import gzip
import json
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch, Mock
from pathlib import Path

import numpy as np
import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "vendor/MFM-DNARNAFM")]
from integrate_wgs_rna.direct_compact_writer import DirectCompactChrWriter, FOURBIT_CODE_TO_ASCII
from integrate_wgs_rna.window_index import WindowIndexConfig, run_window_index_generation, reusable_metadata
from integrate_wgs_rna.window_policy import persist_window_annotations, validate_window_annotations
from pre_training.preprocess_raw_readset.config import ExportConfig
from pre_training.preprocess_raw_readset.shared_hap_writer import prepare_hap_context
from pre_training.preprocess_raw_readset.windowing import sequence_annotations
from rna_window_coverage import retained_windows


class AnnotatedWindowTest(unittest.TestCase):
    def test_n_annotations_preserve_values_without_selection_flags(self):
        for sequence, count, fraction, longest in (
            ("A" * 100, 0, 0.0, 0),
            ("N" + "A" * 99, 1, 0.01, 1),
            ("NN" + "A" * 98, 2, 0.02, 2),
            ("N" + "A" * 24, 1, 0.04, 1),
            ("N" * 100, 100, 1.0, 100),
        ):
            with self.subTest(length=len(sequence), n=sequence.count("N")):
                observed = sequence_annotations(sequence)
                self.assertEqual((observed["n_count"], observed["n_fraction"], observed["n_max_run"]),
                                 (count, fraction, longest))
                self.assertFalse({"has_n", "n_gt_1pct", "training_n5_eligible"}.intersection(observed))

    def test_legacy_exporter_accepts_default_internal_n(self):
        from pre_training.preprocess_raw_readset.exporter import export_readset
        with tempfile.TemporaryDirectory() as raw:
            base=Path(raw)
            sequence="ACGNACGT"+"N"*8+"ACGT"
            fasta=base/"chr1.fa"
            fasta.write_text(">chr1\n"+sequence+"\n")
            pysam.faidx(str(fasta))
            chain=base/"chr1.chain"
            chain.write_text(f"chain 0 chr1 {len(sequence)} + 0 {len(sequence)} chr1 {len(sequence)} + 0 {len(sequence)} 1\n{len(sequence)}\n")
            bam=base/"empty.bam"
            with pysam.AlignmentFile(str(bam),"wb",header={"HD":{"VN":"1.6","SO":"coordinate"},"SQ":[{"SN":"chr1","LN":len(sequence)}]}):
                pass
            pysam.index(str(bam))
            qc=export_readset(sample="S1",contig="chr1",individual_fasta=fasta,bam_path=bam,chain_path=chain,
                              output_dir=base/"export",config=ExportConfig(window_size=8,offset_unit=1))
            self.assertEqual(qc["status"],"pass")
            self.assertEqual(qc["window_count"],2)
            legacy=export_readset(sample="S1",contig="chr1",individual_fasta=fasta,bam_path=bam,chain_path=chain,
                                  output_dir=base/"legacy",config=ExportConfig(window_size=8,offset_unit=1,max_n_ratio=0.0))
            self.assertEqual(legacy["status"],"pass")
            self.assertEqual(legacy["window_count"],1)

    def test_new_default_roundtrip_and_final_mapping(self):
        with tempfile.TemporaryDirectory() as raw:
            base = Path(raw)
            sample = base / "S1"
            sequence = "NNN" + "ACGT" * 25 + "A" * 45 + "N" + "A" * 54 + "N" * 100 + "A" * 12 + "N" + "A" * 12 + "NNNN"
            for hap in ("HAP1", "HAP2"):
                directory = sample / "individual_fasta" / hap
                directory.mkdir(parents=True)
                fasta = directory / "chr1.fa"
                fasta.write_text(">chr1\n" + sequence + "\n")
                pysam.faidx(str(fasta))
                (directory / "chr1.chain").write_text(f"chain 0 chr1 {len(sequence)} + 0 {len(sequence)} chr1 {len(sequence)} + 0 {len(sequence)} 1\n{len(sequence)}\n")
            bed = base / "blacklist.bed"
            bed.write_text("chr1\t3\t103\n")
            cfg = WindowIndexConfig(cohort="TEST", output_root=base, samples=["S1"], chroms=["chr1"],
                                    sequence_variant_policy="snp_only", window_size=100, offset_unit=1, blacklist_bed=bed)
            run_window_index_generation(cfg)
            index = sample / "windows"
            with (index / "windows.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            one = [r for r in rows if r["haplotype"] == "HAP1"]
            self.assertEqual([int(r["raw_start"]) for r in one if r["kept"] == "1"], [3, 103, 303])
            self.assertEqual(one[0]["bl_frac"], "1.000000")
            self.assertEqual(one[1]["n_count"], "1")
            self.assertEqual(one[2]["storage_drop_reason"], "all_n")
            self.assertEqual(one[3]["n_fraction"], "0.04")
            self.assertTrue(all(not {"has_n", "n_gt_1pct", "training_n5_eligible"}.intersection(r) for r in rows))
            # 구캐시의 파생 flag만 제거함. FASTA 재읽기·원값/좌표/선택 변경 없음.
            before = [dict(r) for r in rows]
            legacy = [{**r, "has_n": "1", "n_gt_1pct": "0", "training_n5_eligible": "1"} for r in rows]
            with (index / "windows.tsv").open("w", newline="") as handle:
                out = csv.DictWriter(handle, delimiter="\t", fieldnames=list(legacy[0]))
                out.writeheader()
                out.writerows(legacy)
            old_meta = json.loads((index / "windows.json").read_text())
            old_meta["n_annotation_schema"] = "mfm_window_n_annotations_v2"
            old_meta["n_flag_definitions"] = {"has_n": "legacy", "n_gt_1pct": "legacy"}
            old_meta["training_max_n_fraction"] = 0.05
            (index / "windows.json").write_text(json.dumps(old_meta))
            with patch("integrate_wgs_rna.window_index.read_fasta_sequence", side_effect=AssertionError("FASTA 재읽기 금지")):
                reused = run_window_index_generation(cfg)
                self.assertEqual(reused["skip_count"], 1)
                first_upgrade = (index / "windows.tsv").read_bytes()
                run_window_index_generation(cfg)
                self.assertEqual((index / "windows.tsv").read_bytes(), first_upgrade)
            with (index / "windows.tsv").open() as handle:
                upgraded = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(before, upgraded)
            refreshed_meta = json.loads((index / "windows.json").read_text())
            self.assertNotIn("training_max_n_fraction", refreshed_meta)
            self.assertNotIn("n_flag_definitions", refreshed_meta)
            self.assertEqual(refreshed_meta["n_annotation_schema"], "mfm_window_n_annotations_v3")
            layout = sample / "compact/layout"
            target = layout / "chr1"
            writer = DirectCompactChrWriter(target, final_target=target, sample="S1", contig="chr1",
                                            hap_window_counts={"HAP1": 3, "HAP2": 3}, window_size=100,
                                            source_alignment=base / "source.bam", reference_fasta=None, sample_metadata={},
                                            attrs={"read_coordinate_mode": "personal_realigned"})
            contexts = {}
            for hap in ("HAP1", "HAP2"):
                directory = sample / "individual_fasta" / hap
                context = prepare_hap_context(sample="S1", label=hap, contig="chr1", individual_fasta=directory / "chr1.fa",
                                              chain_path=directory / "chr1.chain", config=ExportConfig(window_size=100, offset_unit=1),
                                              window_index_path=index)
                windows = [row[2] for row in context["eligible_rows"]]
                contexts[hap] = windows
                writer.append_hap_chunk(hap, {"windows": windows, "records": [[] for _ in windows], "attrs": {}, "counts": {}})
            writer.finalize(validation_level="full")
            group = writer.hap_arrays["HAP1"]
            packed = np.asarray(group["bases_packed"][:], dtype=np.uint8)
            codes = np.empty(packed.size * 2, dtype=np.uint8)
            codes[0::2], codes[1::2] = packed >> 4, packed & 15
            decoded = FOURBIT_CODE_TO_ASCII[codes].tobytes().decode()
            self.assertEqual(decoded, "".join(w.sequence for w in contexts["HAP1"]))
            self.assertIn("N", decoded)
            self.assertTrue(decoded.endswith("P" * 75))
            rna = retained_windows(index / "windows.tsv", "HAP1", {"chr1": len(sequence)}, ["chr1"])
            self.assertEqual([r["personal_start"] for r in rna], [3, 103, 303])
            record = {"contig": "chr1", "hap1_chain": str(sample / "individual_fasta/HAP1/chr1.chain"),
                      "hap2_chain": str(sample / "individual_fasta/HAP2/chr1.chain")}
            # 기존 compact에 남은 flag도 최종 보존 단계에서 제거함.
            with (target / "windows.tsv").open() as handle:
                old_compact = list(csv.DictReader(handle, delimiter="\t"))
            old_compact = [{**r, "has_n": "1", "n_gt_1pct": "0", "training_n5_eligible": "1"} for r in old_compact]
            with (target / "windows.tsv").open("w", newline="") as handle:
                out = csv.DictWriter(handle, delimiter="\t", fieldnames=list(old_compact[0]))
                out.writeheader()
                out.writerows(old_compact)
            summary = persist_window_annotations(sample, layout, index, [record], {"test": True})
            self.assertEqual(summary["stored_windows"], 6)
            self.assertEqual(validate_window_annotations(sample)["status"], "pass")
            with (target / "windows.tsv").open() as handle:
                final_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(final_rows[1]["n_count"], "1")
            self.assertTrue(all(not {"has_n", "n_gt_1pct", "training_n5_eligible"}.intersection(r) for r in final_rows))
            self.assertEqual(final_rows[2]["n_fraction"], "0.04")
            self.assertEqual(final_rows[2]["actual_length"], "25")
            self.assertEqual(summary["n_annotation_schema"], "mfm_window_n_annotations_v3")
            self.assertNotIn("max_n_fraction", summary["training_selection"])
            # 실제 full_preprocess 진입점이 주석 보존까지 연결되는지 검증함.
            # 무거운 정렬/writer만 대체하고 기존 Zarr·window·chain은 실제 사용함.
            from integrate_wgs_rna.full_preprocess import FullPreprocessConfig, run_full_preprocess
            config = FullPreprocessConfig(
                sample_id="S1", dataset_role="pre_training", read_type="short_read", project_id="TEST",
                hap1_fasta_dir=sample / "individual_fasta/HAP1", hap1_chain_dir=sample / "individual_fasta/HAP1",
                hap2_fasta_dir=sample / "individual_fasta/HAP2", hap2_chain_dir=sample / "individual_fasta/HAP2",
                alignment=base / "source.bam", contigs=["chr1"], output_sample_root=sample,
                window_index_dir=index, layout_name="layout", storage_policy="keep_raw_intermediate",
                window_size=100, offset_unit=1, sequence_variant_policy="snp_only")
            patches = {
                "validate_read_coordinate_mode": Mock(),
                "validate_personal_fasta_reference_contract": Mock(return_value={"test": True}),
                "build_records": Mock(return_value=[record]),
                "inspect_alignment_phase_tags": Mock(return_value={"hp_present_count": 1}),
                "run_direct_compact_jobs": Mock(return_value={"status": "complete"}),
                "summarize_direct_compact_outputs": Mock(return_value={}),
                "write_full_metadata": Mock(),
            }
            with patch.multiple("integrate_wgs_rna.full_preprocess", **patches):
                completed = run_full_preprocess(config)
            self.assertEqual(completed["status"], "complete")
            self.assertEqual(validate_window_annotations(sample)["stored_windows"], 6)
            # 실제 정리 후 원본 FASTA/window 없이도 최종 계약을 검증한다.
            for hap in ("HAP1", "HAP2"):
                shutil.copytree(target / "readset.zarr", target / f"readset.{hap}.zarr")
            (layout / "manifest.json").write_text(json.dumps({"records": {"chr1": {"contig": "chr1", "status": "pass"}}}))
            (sample / "qc").mkdir(exist_ok=True)
            shutil.rmtree(sample / "individual_fasta")
            shutil.rmtree(index)
            shutil.rmtree(sample / "work", ignore_errors=True)
            from handoff.final_check import validate
            final = validate(sample, "S1", "wgs")
            self.assertEqual(final["status"], "pass")
            self.assertEqual(final["window_annotations_and_mapping"]["stored_windows"], 6)
            mapping = sample / summary["mappings"][0]["path"]
            with mapping.open("ab") as handle:
                handle.write(b"corrupt")
            with self.assertRaisesRegex(ValueError, "SHA256"):
                validate_window_annotations(sample)
            with self.assertRaisesRegex(ValueError, "mapping"):
                validate(sample, "S1", "wgs")


if __name__ == "__main__":
    unittest.main()
