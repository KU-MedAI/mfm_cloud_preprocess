from __future__ import annotations

import hashlib
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pysam
import zarr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from integrate_wgs_rna import direct_compact_writer as writer
from integrate_wgs_rna import stage70_trace_validator as trace


def _fasta(path: Path, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">chr1\n{sequence}\n", encoding="ascii")
    pysam.faidx(str(path))


def _cram(path: Path, target: Path, sequence: str) -> None:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": len(sequence), "M5": hashlib.md5(sequence.encode()).hexdigest()}]}
    read = pysam.AlignedSegment()
    read.query_name, read.query_sequence, read.query_qualities = "r1", "ACGT", pysam.qualitystring_to_array("IIII")
    read.flag, read.reference_id, read.reference_start, read.mapping_quality, read.cigar = 0, 0, 0, 60, ((0, 4),)
    with pysam.AlignmentFile(str(path), "wc", header=header, reference_filename=str(target)) as out:
        out.write(read)
    pysam.index(str(path))


class Stage70TraceValidatorTest(unittest.TestCase):
    def test_symlink_CRAM은원래경로의sidecarindex를사용함(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);ref=root/'ref.fa';_fasta(ref,'ACGT'*4)
            actual=root/'cache/data.cram';actual.parent.mkdir();_cram(actual,ref,'ACGT'*4)
            logical=root/'published.cram';logical.symlink_to(actual)
            index=Path(str(actual)+'.crai');shutil.copy2(index,Path(str(logical)+'.crai'));index.unlink()
            metadata={'source_alignment':{'HAP1':str(logical)},'input_provenance':{'window_index':{'path':str(ref)},'haplotypes':{'HAP1':{'alignment':{'path':str(actual)},'target_fasta':{'path':str(ref)},'hap_fasta':{'path':str(ref)}}}}}
            chosen,*_=trace._paths(metadata,'HAP1');self.assertEqual(chosen,logical)
            with pysam.AlignmentFile(str(chosen),'rb',reference_filename=str(ref)) as handle:self.assertEqual(len(list(handle.fetch('chr1'))),1)
            other=root/'other.cram';_cram(other,ref,'ACGT'*4);metadata['source_alignment']['HAP1']=str(other)
            with self.assertRaisesRegex(ValueError,'provenance mismatch'):trace._paths(metadata,'HAP1')

    def test_sqlite_multiset_batches_preserve_multiplicity(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            bag = trace._Bag(Path(raw) / "bag.sqlite3", batch_size=2)
            try:
                self.assertEqual(bag.db.execute("PRAGMA journal_mode").fetchone()[0], "off")
                self.assertEqual(bag.db.execute("PRAGMA synchronous").fetchone()[0], 0)
                self.assertEqual(bag.db.execute("PRAGMA locking_mode").fetchone()[0], "exclusive")
                digest_a = hashlib.sha256(b"a").digest()
                digest_b = hashlib.sha256(b"b").digest()
                for side in ("src", "out"):
                    bag.add(side, "read", digest_a)
                    bag.add(side, "read", digest_a)
                    bag.add(side, "read", digest_b)
                self.assertTrue(bag.result(2)["pass"])
                bag.add("src", "read", digest_b)
                mismatch = bag.result(2)
                self.assertFalse(mismatch["pass"])
                self.assertEqual(mismatch["mismatches"]["read"]["source_only"], 1)
            finally:
                bag.close()

    def test_root_forwards_requested_scratch_directory(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            compact_root = root / "compact"
            for contig in ("chr2", "chr1"):
                path = compact_root / contig
                path.mkdir(parents=True)
                (path / "contig_metadata.json").write_text(json.dumps({
                    "sample_id": "TEST", "contig": contig,
                    "read_coordinate_mode": "personal_realigned", "readset_layout": "hap_split",
                    "input_provenance": {"fingerprint_sha256": f"fingerprint-{contig}"},
                }), encoding="utf-8")
            observed: list[tuple[str, Path]] = []

            def fake_validate(compact, metadata, hap, scratch, max_examples):
                observed.append((hap, Path(scratch)))
                return {"status": "pass"}

            original = trace._safe_validate_hap
            try:
                trace._safe_validate_hap = fake_validate
                result = trace.validate_stage70_trace_root(
                    compact_root, ["chr2", "chr1"], root / "report",
                    scratch_dir=root / "scratch", workers=1,
                )
            finally:
                trace._safe_validate_hap = original
            self.assertEqual(result["summary"]["status"], "pass")
            self.assertEqual(list(result["results"]), ["chr2", "chr1"])
            self.assertEqual(observed, [
                ("HAP1", root / "scratch"), ("HAP2", root / "scratch"),
                ("HAP1", root / "scratch"), ("HAP2", root / "scratch"),
            ])
            self.assertTrue((root / "report" / "chr1.trace_qc.json").is_file())

    def test_source_zarr_trace_passes_then_detects_read_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            windows = root / "windows"
            windows.mkdir()
            (windows / "windows.tsv").write_text(
                "contig\thaplotype\traw_start\tactual_length\tkept\n"
                "chr1\tHAP1\t0\t8\t1\nchr1\tHAP1\t8\t8\t1\n"
                "chr1\tHAP2\t0\t8\t1\nchr1\tHAP2\t8\t8\t1\n",
                encoding="utf-8",
            )
            job: dict[str, object] = {
                "contig": "chr1", "read_coordinate_mode": "personal_realigned",
                "validation_level": "full", "compact_compressor": "zstd", "compact_compression_level": 5,
                "config": {
                    "window_size": 8, "offset_unit": 1, "curriculum_unit": 8, "seed": 42,
                    "min_mapping_quality": 30, "exclude_unmapped": True, "exclude_secondary": True,
                    "exclude_supplementary": False, "exclude_qcfail": True, "exclude_duplicate": False,
                    "require_proper_pair": False, "max_n_ratio": 0.0, "min_coordinate_coverage": 0.0,
                },
                "window_index_dir": str(windows), "max_windows": 1,
            }
            for hap in (1, 2):
                fasta = root / "individual" / f"HAP{hap}" / "chr1.fa"
                target = root / "realigned" / f"HAP{hap}" / "target.fa"
                cram = root / "realigned" / f"HAP{hap}" / f"HAP{hap}.cram"
                _fasta(fasta, "ACGTACGT")
                _fasta(target, "ACGTACGT")
                _cram(cram, target, "ACGTACGT")
                job[f"hap{hap}_fasta"] = str(fasta)
                job[f"hap{hap}_alignment"] = str(cram)
            provenance = writer._personal_input_provenance(job, {})
            self.assertEqual(provenance["max_windows"], 1)
            compact_root = root / "compact"
            compact_root.mkdir()
            for label in ("HAP1", "HAP2"):
                number = label[-1]
                alignment = Path(str(job[f"hap{number}_alignment"]))
                single = root / f".{label}"
                compact = writer.DirectCompactChrWriter(
                    single, final_target=single, sample="TEST", contig="chr1", hap_window_counts={label: 1},
                    window_size=8, source_alignment=alignment, reference_fasta=alignment.parent / "target.fa",
                    sample_metadata={}, attrs={"read_coordinate_mode": "personal_realigned"},
                )
                compact.append_hap_chunk(label, {
                    "windows": [SimpleNamespace(sequence="ACGTACGT", raw_start=0, actual_length=8)],
                    "records": [[{
                        "read_key": ("r1", 0), "shared_read_key": ("r1", 0, 0, 4, "4M", False, False),
                        "read_sequence": "ACGT", "read_haplotype": 0, "read_phase_set": 0,
                        "query_start": 0, "query_end": 4, "read_start": 0, "read_end": 4,
                        "read_strand": 1, "reference_start": 0, "reference_end": 4,
                        "cigar_ops_packed": [4 << 4], "segment_type": 0,
                    }]], "attrs": {}, "counts": {},
                })
                compact.finalize(validation_level="full")
                if label == "HAP1":
                    shutil.move(str(single / "hap_sequence.zarr"), str(compact_root / "hap_sequence.zarr"))
                    shutil.move(str(single / "depth.zarr"), str(compact_root / "depth.zarr"))
                else:
                    shutil.move(str(single / "hap_sequence.zarr" / label), str(compact_root / "hap_sequence.zarr" / label))
                    shutil.move(str(single / "depth.zarr" / label), str(compact_root / "depth.zarr" / label))
                shutil.move(str(single / "readset.zarr"), str(compact_root / f"readset.{label}.zarr"))
                shutil.rmtree(single)
            (compact_root / "contig_metadata.json").write_text(json.dumps({
                "sample_id": "TEST", "contig": "chr1", "read_coordinate_mode": "personal_realigned",
                "readset_layout": "hap_split", "input_provenance": provenance,
            }), encoding="utf-8")
            passed = trace.validate_stage70_trace_contig(compact_root, scratch_dir=root, max_examples=2, workers=2)
            self.assertEqual(passed["status"], "pass", json.dumps(passed, indent=2))
            # 독립 trace도 SEQ/CIGAR/read metadata/segment/depth 변조를 모두 탐지해야 한다.
            readset = zarr.open_group(str(compact_root / "readset.HAP1.zarr"), mode="a")
            depth = zarr.open_group(str(compact_root / "depth.zarr"), mode="a")
            for name, array in (
                ("read_metadata", readset["reads"]["ref_start"]),
                ("sequence", readset["reads"]["bases_packed"]),
                ("cigar", readset["cigar"]["ops_packed"]),
                ("segment", readset["segments"]["HAP1"]["query_start"]),
                ("depth", depth["HAP1"]["depth"]),
            ):
                with self.subTest(corruption=name):
                    original = int(array[0])
                    array[0] = original + 1
                    try:
                        failed = trace.validate_stage70_trace_contig(compact_root, scratch_dir=root, workers=2)
                        self.assertEqual(failed["status"], "fail", name)
                        hap_qc = failed["by_haplotype"]["HAP1"]
                        if name == "depth":
                            self.assertEqual(hap_qc["projection"]["depth_mismatch_bases"], 1)
                        else:
                            self.assertFalse(hap_qc["record_multiset"]["pass"])
                    finally:
                        array[0] = original
            self.assertEqual(trace.validate_stage70_trace_contig(compact_root, scratch_dir=root)["status"], "pass")


if __name__ == "__main__":
    unittest.main()
