from __future__ import annotations

import hashlib
import shlex
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np
import pysam
import zarr


ROOT = Path(__file__).resolve().parents[1]


class RnaWindowPackageTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name)
        self.sample_root = self.base / "S1"
        self.layout = self.sample_root / "compact" / "layout"
        self.windows = self.sample_root / "windows" / "windows.tsv"
        self.windows.parent.mkdir(parents=True)
        self.windows.write_text(
            "contig\thaplotype\twindow_index\traw_start\tactual_length\tpersonal_start\tpersonal_end\tref_start\tref_end\tkept\n"
            "chr1\tHAP1\t0\t0\t100\t0\t100\t0\t100\t1\n"
            "chr1\tHAP1\t1\t0\t100\t0\t100\t0\t100\t0\n"
            "chr1\tHAP2\t0\t0\t100\t0\t100\t0\t100\t1\n"
            "chr1\tHAP2\t1\t0\t100\t0\t100\t0\t100\t0\n"
        )
        self.fai = self.base / "target.fa.fai"
        self.fai.write_text("chr1\t100\t6\t100\t101\n")
        self.bam = self.base / "rna.bam"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100}]}
        with pysam.AlignmentFile(self.bam, "wb", header=header) as out:
            read = pysam.AlignedSegment()
            read.query_name = "r1"
            read.query_sequence = "A" * 50
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 10
            read.mapping_quality = 255
            read.cigar = ((0, 50),)
            out.write(read)
        pysam.index(str(self.bam))

    def tearDown(self) -> None:
        self.temp.cleanup()

    def make_coverage(self, hap: str, directory: Path) -> None:
        # 하나의 Stage 18 실행이 전장과 WGS-window Zarr를 함께 만든다.
        subprocess.run(
            [
                sys.executable, str(ROOT / "18_coverage_zarr.py"),
                "--bam", str(self.bam), "--fai", str(self.fai),
                "--out", str(directory / "coverage.zarr"),
                "--windows-tsv", str(self.windows),
                "--window-out", str(directory / "coverage_windows.zarr"),
                "--contigs", "chr1", "--sample", "S1", "--hap", hap,
                "--samtools", "samtools",
            ],
            check=True, capture_output=True, text=True,
        )

    def build_package(self) -> tuple[Path, Path]:
        (self.sample_root / "metadata").mkdir(parents=True)
        (self.sample_root / "metadata" / "sample.json").write_text(json.dumps({
            "sample_id": "S1", "project_id": "test", "cohort": "test", "layout_name": "layout",
            "hap_read_mode": "all_to_both", "rna_personal": {"rna_personal_root": "x"},
        }))
        qc = self.sample_root / "qc"
        qc.mkdir()
        (qc / "wgs_preflight.json").write_text(json.dumps({
            "schema": "mfm_preflight_qc_v2", "sample": "S1", "status": "pass",
        }))
        (qc / "rna_preflight.json").write_text(json.dumps({
            "schema": "mfm_rna_preflight_contract_v2", "sample": "S1", "status": "pass",
        }))
        logs = self.sample_root / "logs" / "stage_logs"
        logs.mkdir(parents=True)
        provenance = []
        for branch, stages in (("WGS", range(1, 13)), ("RNA", range(13, 25))):
            branch_dir = logs / branch
            branch_dir.mkdir()
            for stage in stages:
                label = f"{stage:02d}"
                log_path = branch_dir / f"{label}.log"
                log_path.write_text("ok\n")
                provenance.append({
                    "schema": "mfm_execution_provenance_v1", "event": "finish",
                    "sample_id": "S1", "branch": branch, "stage": label,
                    "stage_name": "fixture", "status": "success", "exit_code": 0,
                    "log_path": str(log_path),
                })
        (self.sample_root / "metadata" / "execution_provenance.jsonl").write_text(
            "".join(json.dumps(row) + "\n" for row in provenance)
        )

        # final validator가 cleanup에 남길 개인 HAP FASTA/FAI/chain을 직접 검사한다.
        reference = self.base / "reference.fa"
        reference.write_text(">chr1\n" + "A" * 100 + "\n")
        Path(str(reference) + ".fai").write_text("chr1\t100\t6\t100\t101\n")
        assets = self.sample_root / "individual_fasta"
        assets.mkdir()
        (assets / "personal_fasta_qc.json").write_text(json.dumps({
            "reference_fasta": str(reference), "sequence_variant_policy": "snv_indel",
        }))
        for hap in ("HAP1", "HAP2"):
            directory = assets / hap
            directory.mkdir()
            (directory / "chr1.fa").write_text(">chr1\n" + "A" * 100 + "\n")
            (directory / "chr1.fa.fai").write_text("chr1\t100\t6\t100\t101\n")
            (directory / "chr1.chain").write_text(
                "chain 0 chr1 100 + 0 100 chr1 100 + 0 100 1\n100\n"
            )

        contig = self.layout / "chr1"
        contig.mkdir(parents=True)
        for name in ("contig_metadata.json", "windows.tsv"):
            (contig / name).write_text("{}\n")
        (contig / "windows.tsv").write_text(
            "contig\thaplotype\twindow_index\traw_start\tactual_length\n"
            "chr1\tHAP1\t0\t0\t100\nchr1\tHAP2\t0\t0\t100\n"
        )
        for name in ("hap_sequence.zarr", "depth.zarr", "readset.HAP1.zarr", "readset.HAP2.zarr"):
            (contig / name).mkdir()
        (self.layout / "manifest.json").write_text(json.dumps({
            "status": "complete", "validation_level": "full",
            "records": {"chr1": {"contig": "chr1", "status": "pass"}},
        }))

        rna = self.sample_root / "rna_personal" / "panel"
        rna.mkdir(parents=True)
        hap_variants = rna / "hap_variants.tsv"
        hap_variants.write_text(
            "contig\tstatus\tkind\th1_start\th1_end\th2_start\th2_end\tgenes\n"
            "chr1\tvalidated\tSNV\t10\t11\t10\t11\tG1\n"
        )
        hv_sha = hashlib.sha256(hap_variants.read_bytes()).hexdigest()
        (rna / "hap_variants_qc.json").write_text(json.dumps({
            "schema": "mfm_hap_variants_qc_v1", "output_sha256": hv_sha,
        }))
        for name in (
            "allelic_counts.tsv", "allelic_counts.split_half.tsv",
            "allelic_counts.both_haps.tsv", "allelic_counts.mmap_adjusted.tsv",
            "ase_training_eligibility.tsv", "rna_realign_summary.json",
        ):
            (rna / name).write_text("x\n")
        count_paths = {
            "counts": rna / "allelic_counts.tsv",
            "split_half": rna / "allelic_counts.split_half.tsv",
            "both_haps": rna / "allelic_counts.both_haps.tsv",
            "mmap_adjusted": rna / "allelic_counts.mmap_adjusted.tsv",
        }
        count_hashes = {name: hashlib.sha256(path.read_bytes()).hexdigest()
                        for name, path in count_paths.items()}
        count_meta = rna / "allelic_counts.meta.json"
        count_meta.write_text(json.dumps({
            "schema": "mfm_allelic_counts_meta_v5",
            "script_revision": "2026-10-02_validated_hap_variant_evidence",
            "hap_variants_sha256": hv_sha,
            "output_sha256": count_hashes,
        }))
        (rna / "ase_qc.json").write_text(json.dumps({
            "schema": "mfm_ase_qc_v3",
            "script_revision": "2026-10-02_validated_hap_variant_evidence_cnv_hap_consistency",
            "inputs": {
                "hap_variants_sha256": hv_sha,
                "counts_sha256": count_hashes["counts"],
                "split_half_sha256": count_hashes["split_half"],
                "both_haps_sha256": count_hashes["both_haps"],
                "mmap_adjusted_sha256": count_hashes["mmap_adjusted"],
                "meta_sha256": hashlib.sha256(count_meta.read_bytes()).hexdigest(),
            },
            "sample_qc": {"provenance_ok": True, "duplicates_excluded": True},
        }))
        for hap in ("HAP1", "HAP2"):
            directory = rna / hap
            directory.mkdir()
            for name in ("annotation.gtf", "exons.bygene.tsv", "gtf_lift_qc.json", "rna_realign.json"):
                (directory / name).write_text("x\n")
            self.make_coverage(hap, directory)

        trace = qc / "stage12_trace" / "all_to_both" / "trace_manifest.json"
        trace.parent.mkdir(parents=True)
        time.sleep(0.01)
        trace.write_text(json.dumps({
            "status": "pass", "compact_root": str(self.layout.resolve()), "contigs": {"chr1": "pass"},
        }))
        output = qc / "final_package_validation.json"
        complete = self.sample_root / "metadata" / "completeness.tsv"
        return output, complete

    def validate(self, output: Path, complete: Path, signal=False) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable, str(ROOT / "scripts/validate_sample_package.py"),
                "--sample", "S1", "--root", str(self.sample_root), "--mode", "rna",
                "--layout", "layout", "--chroms", "chr1", "--hap-read-mode", "all_to_both",
                "--rna-label", "panel", "--trace-policy", "required",
                "--rna-payload", "coverage_only" if signal else "full",
                "--output", str(output), "--completeness-tsv", str(complete),
            ], capture_output=True, text=True,
        )

    def test_signal_only_full_validator_without_ase_and_depth_tamper(self):
        from rna_signal import stamp_signal
        from scripts.cleanup_guard import validate_predelete
        output, complete = self.build_package()
        directory=self.sample_root/'rna_personal/panel'
        for path in directory.iterdir():
            if path.is_file() and path.name.startswith(('hap_variants','allelic_counts','ase_','gene_rd')):path.unlink()
        preflight_path=self.sample_root/'qc/rna_preflight.json'
        preflight=json.loads(preflight_path.read_text());content={k:v for k,v in preflight.items() if k not in ('contract_sha256','status')}
        digest=hashlib.sha256(json.dumps(content,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        preflight['contract_sha256']=digest;preflight_path.write_text(json.dumps(preflight))
        fastq=dict(schema='mfm_rna_fastq_qc_v2',rna_preflight_contract_sha256=digest)
        haps=[]
        for hap in ('HAP1','HAP2'):
            qc=dict(schema='mfm_v2_rna_realign_v6',sample='S1',hap=hap,phase_mode='panel',fastq_reconstruction=fastq)
            (directory/hap/'rna_realign.json').write_text(json.dumps(qc));haps.append(qc)
        summary=dict(schema='mfm_rna_qc_v2',sample='S1',phase_mode='panel',haplotypes=haps,fastq_reconstruction=fastq)
        for name in ('rna_qc.json','rna_realign_summary.json'):(directory/name).write_text(json.dumps(summary))
        stamp_signal(self.sample_root)
        # 실제 공개 셸이 RNA 모드를 최종 validator까지 전달하는지 확인함.
        entry=self.base/'entry';entry.mkdir()
        wrapper=ROOT/('25_validate_sample_package.sh' if (ROOT/'09_personal_fasta.sh').exists() else '23_validate_sample_package.sh')
        (entry/'validate.sh').write_text(wrapper.read_text())
        (entry/'scripts').symlink_to(ROOT/'scripts',target_is_directory=True)
        config="PYTHON="+shlex.quote(sys.executable)+"\nCOMPACT_LAYOUT_NAME=layout\nCHROMS=chr1\nHAP_READ_MODE=all_to_both\nSTAGE12_TRACE_POLICY=required\nMFM_RNA_OUTPUT_MODE=coverage_only\n"
        config+="apply_sample_payload() { :; }\nchroms_space() { echo chr1; }\nsample_rna_evidence_label() { echo panel; }\n"
        config+="sample_final_root() { echo "+shlex.quote(str(self.sample_root.parent))+"; }\n"
        config+="sample_work_root() { echo "+shlex.quote(str(self.base/'no_work'))+"; }\n"
        (entry/'config_v2.sh').write_text(config)
        first=subprocess.run(['bash',str(entry/'validate.sh'),'S1','rna'],capture_output=True,text=True)
        self.assertEqual(first.returncode,0,first.stderr+first.stdout)
        self.assertEqual(validate_predelete(self.sample_root,'S1','rna')['rna_payload'],'coverage_only')
        group=zarr.open_group(str(directory/'HAP1/coverage.zarr'),mode='a');group['chr1'][0]=99
        second=self.validate(output,complete,signal=True)
        self.assertNotEqual(second.returncode,0)
        self.assertIn('rna_HAP1:genome_coverage',json.loads(output.read_text())['failed_checks'])

    def test_exact_window_contract_and_tamper_detection(self) -> None:
        output, complete = self.build_package()
        first = self.validate(output, complete)
        self.assertEqual(first.returncode, 0, first.stderr + first.stdout)
        report = json.loads(output.read_text())
        self.assertEqual(report["status"], "pass")
        for hap in ("HAP1", "HAP2"):
            store = self.sample_root / "rna_personal" / "panel" / hap / "coverage_windows.zarr"
            group = zarr.open_group(str(store), mode="r")
            self.assertEqual(group.attrs["window_count"], 1)
            self.assertEqual(list(group["window_index"][:]), [0])
            self.assertEqual(group["depth_mean"][:].tolist(), [0.5])
            self.assertEqual(group["depth_zero_fraction"][:].tolist(), [0.5])
            self.assertEqual(group["depth_max"][:].tolist(), [1])
            genome = zarr.open_group(str(store.parent / "coverage.zarr"), mode="r")
            np.testing.assert_array_equal(group["depth"][:], genome["chr1"][:])

        store = self.sample_root / "rna_personal" / "panel" / "HAP2" / "coverage_windows.zarr"
        group = zarr.open_group(str(store), mode="a")
        group["window_key"][0] = "chr1:HAP2:0:1-100"
        second = self.validate(output, complete)
        self.assertNotEqual(second.returncode, 0)
        report = json.loads(output.read_text())
        self.assertIn("rna_HAP2:window_contract", report["failed_checks"])

    def test_legacy_coordinates_produce_identical_depth(self) -> None:
        directory = self.base / "modern"
        self.make_coverage("HAP1", directory)
        modern = zarr.open_group(str(directory / "coverage_windows.zarr"), mode="r")
        expected = {name: modern[name][:] for name in modern.array_keys()}
        self.windows.write_text(
            "contig\thaplotype\twindow_index\traw_start\tactual_length\tkept\n"
            "chr1\tHAP1\t0\t0\t100\t1\nchr1\tHAP2\t0\t0\t100\t1\n"
        )
        self.make_coverage("HAP1", self.base / "legacy")
        legacy = zarr.open_group(str(self.base / "legacy/coverage_windows.zarr"), mode="r")
        for name, values in expected.items():
            if name in {"ref_start", "ref_end"}:
                # 구형 입력에 REF 투영이 없다. 개인 좌표를 REF로 복사하지 않는다.
                self.assertEqual(legacy[name][:].tolist(), [-1])
                continue
            np.testing.assert_array_equal(values, legacy[name][:], err_msg=name)
        self.assertEqual(int(legacy["depth"][:].sum()), 50)

    def test_coordinate_array_tampering_is_rejected(self) -> None:
        output, complete = self.build_package()
        group = zarr.open_group(str(self.sample_root / "rna_personal/panel/HAP1/coverage_windows.zarr"), mode="a")
        group["personal_start"][0] = 1
        result = self.validate(output, complete)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("rna_HAP1:window_contract", json.loads(output.read_text())["failed_checks"])

    def test_hap_variant_hash_tampering_is_rejected(self) -> None:
        output, complete = self.build_package()
        path = self.sample_root / "rna_personal/panel/hap_variants.tsv"
        path.write_text(path.read_text().replace("\tG1\n", "\tG2\n"))
        result = self.validate(output, complete)
        self.assertNotEqual(result.returncode, 0)
        failed = json.loads(output.read_text())["failed_checks"]
        self.assertIn("rna_hap_variant_qc_contract", failed)
        self.assertIn("rna_allelic_v5_contract", failed)

    def test_ase_input_hash_mismatch_is_rejected(self) -> None:
        output, complete = self.build_package()
        path = self.sample_root / "rna_personal/panel/ase_qc.json"
        payload = json.loads(path.read_text())
        payload["inputs"]["hap_variants_sha256"] = "0" * 64
        path.write_text(json.dumps(payload))
        result = self.validate(output, complete)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "rna_ase_validated_evidence_contract",
            json.loads(output.read_text())["failed_checks"],
        )

    def test_count_content_tampering_is_rejected(self) -> None:
        output, complete = self.build_package()
        path = self.sample_root / "rna_personal/panel/allelic_counts.tsv"
        path.write_text("tampered\n")
        result = self.validate(output, complete)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "rna_ase_validated_evidence_contract",
            json.loads(output.read_text())["failed_checks"],
        )

    def test_coverage_content_tampering_is_rejected(self) -> None:
        output, complete = self.build_package()
        group = zarr.open_group(
            str(self.sample_root / "rna_personal/panel/HAP1/coverage_windows.zarr"), mode="a"
        )
        group["depth"][10] = int(group["depth"][10]) + 1
        result = self.validate(output, complete)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("rna_HAP1:window_contract", json.loads(output.read_text())["failed_checks"])

    def test_genome_coverage_content_tampering_is_rejected(self) -> None:
        output, complete = self.build_package()
        group = zarr.open_group(
            str(self.sample_root / "rna_personal/panel/HAP1/coverage.zarr"), mode="a"
        )
        group["chr1"][10] = int(group["chr1"][10]) + 1
        result = self.validate(output, complete)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("rna_HAP1:genome_coverage", json.loads(output.read_text())["failed_checks"])

    def test_missing_provenance_stage_is_rejected(self) -> None:
        output, complete = self.build_package()
        path = self.sample_root / "metadata/execution_provenance.jsonl"
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        path.write_text("".join(
            json.dumps(row) + "\n"
            for row in rows if not (row["branch"] == "RNA" and row["stage"] == "24")
        ))
        result = self.validate(output, complete)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("execution_provenance", json.loads(output.read_text())["failed_checks"])

    def test_personal_fasta_tampering_after_manifest_is_rejected(self) -> None:
        output, complete = self.build_package()
        first = self.validate(output, complete)
        self.assertEqual(first.returncode, 0, first.stderr + first.stdout)
        fasta = self.sample_root / "individual_fasta/HAP1/chr1.fa"
        fasta.write_text(fasta.read_text().replace("A", "C", 1))
        second = self.validate(output, complete)
        self.assertNotEqual(second.returncode, 0)
        self.assertIn("personal_hap_assets", json.loads(output.read_text())["failed_checks"])

    def test_preflight_content_tampering_is_rejected(self) -> None:
        output, complete = self.build_package()
        path = self.sample_root / "qc/rna_preflight.json"
        payload = json.loads(path.read_text())
        payload["sample"] = "OTHER"
        path.write_text(json.dumps(payload))
        result = self.validate(output, complete)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("rna_preflight_qc", json.loads(output.read_text())["failed_checks"])


if __name__ == "__main__":
    unittest.main()
