"""RNA signal 전용의 입력·완료·삭제 보호 계약을 검증함."""
import csv
import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import zarr

from handoff import planner
from handoff.resume import sample_downloads
from handoff.sequence_mode import selected, needs_wgs, validate_selection, payload_matches
from rna_signal import stamp_signal, validate_signal
from scripts.cleanup_guard import validate_predelete

ROOT = Path(__file__).resolve().parents[1]


class SignalTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name) / 'S1'
        for part in ('metadata', 'qc', 'individual_fasta', 'rna_personal/public/HAP1', 'rna_personal/public/HAP2'):
            (self.root / part).mkdir(parents=True, exist_ok=True)
        self.write('metadata/sample.json', dict(sample_id='S1', wgs_payload='sequence_only'))
        preflight = dict(schema='mfm_rna_preflight_contract_v2', sample='S1', rna={'size': 1})
        digest = hashlib.sha256(json.dumps(preflight, ensure_ascii=False, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        preflight.update(contract_sha256=digest, status='pass')
        self.write('qc/rna_preflight.json', preflight)
        fastq = dict(schema='mfm_rna_fastq_qc_v2', rna_preflight_contract_sha256=digest)
        self.write('individual_fasta/permanent_assets.json', dict(schema='fixture'))
        haps = []
        for hap in ('HAP1', 'HAP2'):
            qc = dict(schema='mfm_v2_rna_realign_v6', sample='S1', hap=hap, phase_mode='public', fastq_reconstruction=fastq)
            haps.append(qc)
            self.write(f'rna_personal/public/{hap}/rna_realign.json', qc)
            for name in ('annotation.gtf', 'exons.bygene.tsv', 'gtf_lift_qc.json'):
                (self.root / 'rna_personal/public' / hap / name).write_text('fixture\n')
            for name in ('coverage.zarr', 'coverage_windows.zarr'):
                g = zarr.open_group(str(self.root / 'rna_personal/public' / hap / name), mode='w')
                g.attrs.update(sample='S1', hap=hap, source_bam_signature={'size': 1}, source_fai_signature={'size': 1})
                g.create_dataset('depth', data=[1, 2, 3], dtype='u4')
        summary = dict(schema='mfm_rna_qc_v2', sample='S1', phase_mode='public', haplotypes=haps, fastq_reconstruction=fastq)
        for name in ('rna_qc.json', 'rna_realign_summary.json'):
            self.write('rna_personal/public/' + name, summary)

    def tearDown(self):
        self.temp.cleanup()

    def write(self, relative, value):
        (self.root / relative).write_text(json.dumps(value))

    def sha(self, relative):
        return hashlib.sha256((self.root / relative).read_bytes()).hexdigest()

    def report(self):
        self.write('qc/final_package_validation.json', dict(
            schema='mfm_final_package_validation_v1', contract_revision='2026-10-03_integrity_v2',
            sample='S1', mode='rna', status='pass', failed_checks=[], rna_payload='coverage_only',
            personal_assets_manifest_sha256=self.sha('individual_fasta/permanent_assets.json'),
            rna_signal_manifest_sha256=self.sha('rna_personal/public/rna_signal_manifest.json')))

    def test_signal_without_ase_passes_guard_and_preserves_identity(self):
        stamp_signal(self.root)
        self.assertEqual(validate_signal(self.root, 'S1', 'public')['rna_payload'], 'coverage_only')
        self.report()
        self.assertEqual(validate_predelete(self.root, 'S1', 'rna')['rna_payload'], 'coverage_only')
        meta = json.loads((self.root / 'metadata/sample.json').read_text())
        self.assertFalse(meta['ase_available'])
        self.assertFalse(meta['cnv_available'])

    def test_pileup_payload_keeps_rna_scope_and_guard(self):
        self.write('metadata/sample.json',dict(sample_id='S1',wgs_payload='pileup_feature'))
        stamp_signal(self.root)
        meta=json.loads((self.root/'metadata/sample.json').read_text())
        self.assertEqual(meta['completion_scope'],'dna_pileup_feature_rna_coverage_only')
        self.assertEqual(meta['wgs_payload'],'pileup_feature')
        self.report()
        self.assertEqual(validate_predelete(self.root,'S1','rna')['rna_payload'],'coverage_only')

    def test_preflight_formatting_preserves_contract_but_content_change_fails(self):
        path = self.root / 'qc/rna_preflight.json'
        data = json.loads(path.read_text())
        path.write_text(json.dumps(data, indent=2, sort_keys=True))
        stamp_signal(self.root)
        data['rna']['size'] = 2
        self.write('qc/rna_preflight.json', data)
        with self.assertRaisesRegex(ValueError, 'preflight'): stamp_signal(self.root)

    def test_post_cleanup_requires_approval_and_preserved_hashes(self):
        stamp_signal(self.root)
        self.report()
        for hap in ('HAP1','HAP2'):
            (self.root/f'rna_personal/public/{hap}/annotation.gtf').unlink()
        with self.assertRaises(ValueError): validate_signal(self.root,'S1','public')
        with self.assertRaises((ValueError,OSError)): validate_signal(self.root,'S1','public',after_cleanup=True)
        (self.root/'metadata/permanent_assets.json').write_bytes((self.root/'individual_fasta/permanent_assets.json').read_bytes())
        self.write('metadata/cleanup_in_progress.json',dict(schema='mfm_delivery_cleanup_v1',sample='S1',mode='rna',dry_run=False,predelete_guard=dict(rna_payload='coverage_only',validation_report_fresh=True)))
        validate_signal(self.root,'S1','public',after_cleanup=True)
        (self.root/'rna_personal/public/HAP1/exons.bygene.tsv').write_text('changed')
        with self.assertRaises(ValueError):validate_signal(self.root,'S1','public',after_cleanup=True)

    def test_missing_hap_qc_blocks_completion(self):
        (self.root / 'rna_personal/public/HAP2/rna_realign.json').unlink()
        with self.assertRaises(ValueError): stamp_signal(self.root)

    def test_qc_mutation_and_source_mutation_reject_old_pass(self):
        stamp_signal(self.root)
        self.report()
        self.write('qc/rna_preflight.json', {'modified': True})
        with self.assertRaises(ValueError): validate_predelete(self.root, 'S1', 'rna')

    def test_newer_track_rejects_cleanup(self):
        stamp_signal(self.root)
        self.report()
        chunk = self.root / 'rna_personal/public/HAP1/coverage.zarr/new_chunk'
        chunk.write_bytes(b'modified')
        newer = (self.root / 'qc/final_package_validation.json').stat().st_mtime_ns + 1_000_000_000
        os.utime(chunk, ns=(newer, newer))
        with self.assertRaisesRegex(ValueError, '오래됨'): validate_predelete(self.root, 'S1', 'rna')

    def test_existing_full_evidence_cannot_be_relabeled(self):
        (self.root / 'rna_personal/public/allelic_counts.tsv').write_text('old full result')
        with self.assertRaises(ValueError): stamp_signal(self.root)

    def test_long_and_short_signal_exclude_wgs_downloads(self):
        short = next(r for r in planner.read_tsv(ROOT / 'manifests/리드별/전체_short_read.tsv') if r['data_group'] == 'wgs_rna_pair')
        long = next(r for r in planner.read_tsv(ROOT / 'manifests/리드별/전체_long_read.tsv') if r['data_group'] == 'wgs_rna_pair')
        c = dict(data_root=self.temp.name, pipeline_env=dict(MFM_SHORT_READ_OUTPUT_MODE='sequence_only', MFM_LONG_READ_OUTPUT_MODE='sequence_only', MFM_RNA_OUTPUT_MODE='coverage_only'))
        resources = [dict(kind=k, path=k, url='fixture') for k in ('wgs_read', 'wgs_alignment', 'wgs_index', 'rna_alignment', 'sample_vcf')]
        for row in (short, long):
            validate_selection(c, [row])
            self.assertTrue(selected(c, row))
            self.assertFalse(needs_wgs(c, row))
            self.assertEqual([r['kind'] for r in sample_downloads(c, row, resources)], ['rna_alignment', 'sample_vcf'])
            self.write('qc/final_package_validation.json', dict(wgs_payload='sequence_only'))
            self.assertFalse(payload_matches(c, row, self.root))
            self.write('qc/final_package_validation.json', dict(wgs_payload='sequence_only', rna_payload='coverage_only'))
            self.assertTrue(payload_matches(c, row, self.root))
        c['pipeline_env']['MFM_SEQUENCE_DEFER_RNA'] = '1'
        with self.assertRaises(ValueError): validate_selection(c, [short])

    def test_real_long_manifest_binds_without_any_bam_or_read_manifest(self):
        from handoff import cli
        row = planner.read_tsv(ROOT / 'manifests/리드별/전체_long_read.tsv')[0]
        c = dict(data_root=self.temp.name, output_root=self.temp.name, state_root=self.temp.name,
                 run_id='fixture', pipeline_env={'MFM_LONG_READ_OUTPUT_MODE': 'sequence_only',
                                                'MFM_RNA_OUTPUT_MODE': 'coverage_only'})
        refs = {k: str(self.root / (k + '.fa')) for k in ('autosomes', 'ncbi', 'broad', 'illumina')}
        with patch('handoff.references.variant_reference', return_value=refs['autosomes']), patch('handoff.bundled_assets.panel_directory', return_value=self.root):
            manifest, bound = cli.bind_manifest(c, row, [], refs)
        self.assertTrue(manifest.is_file())
        self.assertEqual(bound['alignment_path'], '')
        self.assertEqual(bound['read_files_manifest'], '')
        self.assertEqual(bound['wgs_payload'], 'sequence_only')

    def test_verified_local_vcf_needs_no_download_and_changed_file_is_rejected(self):
        from handoff.sequence_mode import uses_bundled_only
        row = planner.read_tsv(ROOT / 'manifests/리드별/전체_long_read.tsv')[0]
        c = dict(data_root=self.temp.name, pipeline_env={'MFM_LONG_READ_OUTPUT_MODE': 'sequence_only', 'MFM_SEQUENCE_DEFER_RNA': '1'})
        p = Path(self.temp.name) / 'sample.vcf.gz'
        p.write_bytes(b'VCF')
        stat = p.stat()
        records = Path(self.temp.name) / '.download_records'
        records.mkdir()
        key = hashlib.sha256(str(p).encode()).hexdigest()
        (records / (key + '.json')).write_text(json.dumps(dict(status='verified', url='fixture', signature=dict(size=stat.st_size, mtime_ns=stat.st_mtime_ns, ctime_ns=stat.st_ctime_ns))))
        resource = dict(path='sample.vcf.gz', url='fixture', shared=False, kind='sample_vcf')
        with patch('handoff.bundled_assets.inventory', return_value=(self.root, {})), patch('handoff.downloads.sample_resources', return_value=[resource]):
            self.assertTrue(uses_bundled_only(c, row, ROOT))
            p.write_bytes(b'changed')
            self.assertFalse(uses_bundled_only(c, row, ROOT))


    def test_queue_requires_post_cleanup_pass_before_counting_complete(self):
        from types import SimpleNamespace
        from handoff.local_queue import Adapter
        row=dict(planner.read_tsv(ROOT/'manifests/리드별/전체_short_read.tsv')[0],sample_id='S1',data_group='wgs',rna_path='')
        adapter=Adapter.__new__(Adapter)
        adapter.config={'pipeline_env':{'CLEAN_INTERMEDIATES_ON_SUCCESS':'1','MFM_SHORT_READ_OUTPUT_MODE':'sequence_only'}}
        adapter.api=SimpleNamespace(sample_root=lambda *a:self.root,processing_complete=lambda *a:True)
        self.write('qc/final_package_validation.json',{'status':'pass','wgs_payload':'sequence_only'})
        self.assertFalse(adapter.verified(row))
        self.assertTrue(adapter.policy_matches(row))
        adapter.config['pipeline_env']['WINDOW_STORAGE_POLICY']='keep_internal_n'
        self.write('metadata/window_policy.json',{'storage_policy':'different'})
        self.assertFalse(adapter.policy_matches(row))
        self.write('metadata/window_policy.json',{'storage_policy':'keep_internal_n'})
        self.assertTrue(adapter.policy_matches(row))
        self.write('metadata/cleanup_completed.json',dict(schema='mfm_delivery_cleanup_v1',sample='S1',mode='wgs',status='complete',dry_run=False,post_cleanup_validation=dict(schema='mfm_delivery_final_validation_v1',sample='S1',mode='wgs',status='pass',errors=[])))
        self.assertTrue(adapter.verified(row))



if __name__ == '__main__':
    unittest.main()
