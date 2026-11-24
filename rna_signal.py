"""ASE/CNV 없는 RNA signal의 명시적 완료·무결성 계약임."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

SCHEMA = 'mfm_rna_signal_manifest_v1'
EXCLUDED = ('hap_variants.tsv', 'hap_variants_qc.json', 'allelic_counts.tsv',
            'allelic_counts.meta.json', 'ase_qc.json', 'ase_training_eligibility.tsv',
            'gene_rd.tsv', 'gene_rd.tsv.meta.json')
REQUIRED = ('rna_realign_summary.json', 'rna_qc.json') + tuple(
    f'{hap}/{name}' for hap in ('HAP1', 'HAP2')
    for name in ('annotation.gtf', 'exons.bygene.tsv', 'gtf_lift_qc.json', 'rna_realign.json'))


def output_mode(environment=None):
    mode = (os.environ if environment is None else environment).get('MFM_RNA_OUTPUT_MODE', 'full')
    if mode not in ('full', 'coverage_only'):
        raise ValueError('MFM_RNA_OUTPUT_MODE는 full/coverage_only여야 함')
    return mode


def _json(path):
    data = json.loads(Path(path).read_text())
    if not isinstance(data, dict):
        raise ValueError(f'RNA JSON object가 아님: {path}')
    return data


def _sha(path):
    state = hashlib.sha256()
    with Path(path).open('rb') as handle:
        while block := handle.read(8 << 20):
            state.update(block)
    return state.hexdigest()


def _evidence(root, sample, label, *, after_cleanup=False):
    import zarr
    directory = root / 'rna_personal' / label
    if root.name != sample or not directory.is_dir() or directory.is_symlink():
        raise ValueError('RNA signal sample/경로 불일치')
    if any((directory / name).exists() for name in EXCLUDED):
        raise ValueError('coverage_only에 ASE/CNV 산출물이 섞임')
    retained = []
    for name in REQUIRED:
        path = directory / name
        if after_cleanup and name.endswith('/annotation.gtf') and not path.exists():
            continue
        if not path.is_file() or path.is_symlink() or not path.stat().st_size:
            raise ValueError(f'RNA signal 필수 파일 누락: {path}')
        retained.append(name)
    summary = _json(directory / 'rna_realign_summary.json')
    if (summary.get('schema') != 'mfm_rna_qc_v2' or summary.get('sample') != sample
            or summary.get('phase_mode') != label
            or summary != _json(directory / 'rna_qc.json')):
        raise ValueError('RNA signal summary 신원/schema 불일치')
    haps = summary.get('haplotypes', [])
    if len(haps) != 2 or {h.get('hap') for h in haps} != {'HAP1', 'HAP2'}:
        raise ValueError('RNA signal 두 HAP QC 필요함')
    fastq = summary.get('fastq_reconstruction', {})
    preflight = _json(root / 'qc/rna_preflight.json')
    # FASTQ는 파일 바이트가 아닌 rna_preflight.py의 정규화 계약 서명을 저장함.
    content = {k: v for k, v in preflight.items() if k not in ('contract_sha256', 'status')}
    digest = hashlib.sha256(json.dumps(content, ensure_ascii=False, sort_keys=True,
                                      separators=(',', ':')).encode()).hexdigest()
    if (preflight.get('schema') != 'mfm_rna_preflight_contract_v2'
            or preflight.get('sample') != sample or preflight.get('status') != 'pass'
            or preflight.get('contract_sha256') != digest
            or fastq.get('schema') != 'mfm_rna_fastq_qc_v2'
            or fastq.get('rna_preflight_contract_sha256') != digest):
        raise ValueError('RNA signal 원본 preflight 서명 불일치')
    for hap in ('HAP1', 'HAP2'):
        qc = _json(directory / hap / 'rna_realign.json')
        if (qc.get('schema') != 'mfm_v2_rna_realign_v6' or qc.get('sample') != sample
                or qc.get('hap') != hap or qc.get('phase_mode') != label
                or qc != next(h for h in haps if h['hap'] == hap)
                or qc.get('fastq_reconstruction') != fastq):
            raise ValueError('RNA signal HAP QC 불일치')
        for store in ('coverage.zarr', 'coverage_windows.zarr'):
            path = directory / hap / store
            group = zarr.open_group(str(path), mode='r')
            if group.attrs.get('sample') != sample or group.attrs.get('hap') != hap:
                raise ValueError('RNA signal Zarr sample/HAP 불일치')
            if not group.attrs.get('source_bam_signature') or not group.attrs.get('source_fai_signature'):
                raise ValueError('RNA signal 입력 서명 누락')
    return {**{name: _sha(directory / name) for name in retained},
            'source_rna_preflight': _sha(root / 'qc/rna_preflight.json')}


def stamp_signal(root):
    root = Path(root)
    meta_path = root / 'metadata/sample.json'
    meta = _json(meta_path)
    labels = [p for p in (root / 'rna_personal').iterdir() if p.is_dir()]
    if len(labels) != 1:
        raise ValueError('RNA signal evidence 디렉터리는 하나여야 함')
    label = labels[0].name
    files = _evidence(root, root.name, label)
    manifest = dict(schema=SCHEMA, sample=root.name, label=label, rna_payload='coverage_only',
                    ase_available=False, cnv_available=False, files=files)
    from handoff.cli import atomic_json
    atomic_json(labels[0] / 'rna_signal_manifest.json', manifest)
    meta.update(rna_payload='coverage_only', rna_processing_status='complete',
                rna_coverage_available=True, ase_available=False, cnv_available=False,
                completion_scope={'sequence_only':'dna_sequence_rna_coverage_only',
                                  'pileup_feature':'dna_pileup_feature_rna_coverage_only',
                                  'full':'wgs_full_rna_coverage_only'}[meta.get('wgs_payload','full')])
    atomic_json(meta_path, meta)


def validate_signal(root, sample, label, *, after_cleanup=False):
    root = Path(root)
    meta = _json(root / 'metadata/sample.json')
    if (meta.get('rna_payload') != 'coverage_only' or meta.get('rna_processing_status') != 'complete'
            or meta.get('rna_coverage_available') is not True
            or meta.get('ase_available') is not False or meta.get('cnv_available') is not False):
        raise ValueError('RNA signal 메타데이터 계약 불일치')
    manifest_path = root / 'rna_personal' / label / 'rna_signal_manifest.json'
    manifest = _json(manifest_path)
    if after_cleanup:
        report = _json(root / 'qc/final_package_validation.json')
        marker = root / 'metadata/cleanup_in_progress.json'
        if not marker.exists(): marker = root / 'metadata/cleanup_completed.json'
        approval = _json(marker)
        guard = approval.get('predelete_guard', {})
        if (report.get('schema') != 'mfm_final_package_validation_v1'
                or report.get('contract_revision') != '2026-10-03_integrity_v2'
                or report.get('sample') != sample or report.get('mode') != 'rna'
                or report.get('status') != 'pass' or report.get('failed_checks')
                or report.get('rna_payload') != 'coverage_only'
                or report.get('rna_signal_manifest_sha256') != _sha(manifest_path)
                or report.get('personal_assets_manifest_sha256') != _sha(root / 'metadata/permanent_assets.json')
                or approval.get('schema') != 'mfm_delivery_cleanup_v1'
                or approval.get('sample') != sample or approval.get('mode') != 'rna'
                or approval.get('dry_run') is not False
                or guard.get('rna_payload') != 'coverage_only'
                or guard.get('validation_report_fresh') is not True):
            raise ValueError('RNA signal 정리 승인/최종검증 서명 불일치')
    observed = _evidence(root, sample, label, after_cleanup=after_cleanup)
    expected = manifest.get('files', {})
    omitted = set(expected) - set(observed)
    if (manifest.get('schema') != SCHEMA or manifest.get('sample') != sample
            or manifest.get('label') != label or manifest.get('rna_payload') != 'coverage_only'
            or (omitted and (not after_cleanup or not omitted <= {'HAP1/annotation.gtf', 'HAP2/annotation.gtf'}))
            or {k: v for k, v in expected.items() if k not in omitted} != observed):
        raise ValueError('RNA signal manifest/schema/파일 SHA 불일치')
    return manifest
