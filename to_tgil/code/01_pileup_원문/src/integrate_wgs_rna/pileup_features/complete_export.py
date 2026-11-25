"""기존 HAP window 서열과 reference 관측의 생성·재개·독립 대조를 연결함."""
from pathlib import Path
import argparse
import hashlib
import json
import multiprocessing
import os
import shutil
import time
import uuid
from collections import Counter

import numpy as np
import pysam
import zarr
from numcodecs import Blosc, blosc
from pre_training.preprocess_raw_readset.config import ExportConfig
from . import MODE
from .schema import LAYOUT, CHANNELS, SCHEMA, ARRAYS, POLICY, FEATURE_SET, resolve_channels
from .complete_compute import ReferenceCache, projected
from .projection import Projection
from .. import sequence_only
from ..cpu_executor import cpu_results
from ..memory_policy import memory_active
from ..full_preprocess import resolve_stage70_filter_policy
from ..pipeline.checkpoints import atomic_json, digest, file_identity
from ..window_contract import read_windows
from ..window_policy import persist_window_annotations, validate_window_annotations

from .legacy_export import (filter_config, reference_contract, hap_reference_contract, _projection,
                            _inputs as base_inputs, _worker_memory_limits)
from .variants import VariantLookup


def _inputs(root, contig, alignment, reference, alignment_reference, config, variants):
    identities = base_inputs(root, contig, alignment, reference, alignment_reference, config)
    index = next((p for p in [Path(str(variants)+'.tbi'), Path(str(variants)+'.csi')] if p.is_file()), None)
    if index is None:
        raise ValueError('BAF sequence VCF index가 필요함')
    return identities + [file_identity(variants), file_identity(index)]


def _contig(root, layout, contig, rows, window_size, alignment, reference, alignment_reference, read_type, min_mapq, variants, variant_sample):
    blosc.set_nthreads(1)
    root, layout = Path(root), Path(layout)
    variants = Path(variants)
    started = time.monotonic()
    print(f'[PILEUP] sample={root.name} contig={contig} phase=start windows={len(rows)}', flush=True)
    alignment, reference, alignment_reference = map(Path, (alignment, reference, alignment_reference))
    config = filter_config(read_type, min_mapq)
    identities = _inputs(root, contig, alignment, reference, alignment_reference, config, variants)
    record = sequence_only._contig_task(root, layout, contig, rows, window_size, MODE)
    directory = layout/contig
    producer_backend = 'reused'
    metadata = json.loads((directory/'contig_metadata.json').read_text())
    provenance = dict(inputs=identities, config=config.to_dict(), policy=POLICY,
                      channels=CHANNELS, schema=SCHEMA, producer=FEATURE_SET, variant_sample=variant_sample)
    fingerprint = digest(provenance)
    if (directory/'pileup.zarr').exists():
        if metadata.get('pileup', {}).get('fingerprint_sha256') != fingerprint:
            raise ValueError('기존 pileup과 source/policy 지문이 다름; 덮어쓰지 않음')
    else:
        for stale in directory.glob('pileup.partial.*'):
            if stale.is_symlink() or not stale.is_dir():
                raise ValueError('안전하게 정리할 수 없는 pileup 임시 경로')
            shutil.rmtree(stale)
        temporary = directory/('pileup.partial.'+uuid.uuid4().hex)
        try:
            codec = Blosc(cname='zstd', clevel=5, shuffle=Blosc.BITSHUFFLE)
            group = zarr.open_group(str(temporary), mode='w')
            group.attrs.update(schema=SCHEMA, wgs_payload=MODE, sample=root.name, contig=contig,
                               channels=list(CHANNELS), policy=POLICY, feature_set=FEATURE_SET, fingerprint_sha256=fingerprint)
            compact_rows = read_windows(directory/'windows.tsv', compact=True)
            with pysam.AlignmentFile(str(alignment), 'r', reference_filename=str(alignment_reference)) as bam:
                contract = reference_contract(bam, reference, alignment_reference, contig)
                hap_reference_contract(root, bam, reference, contig)
                cache = ReferenceCache(bam, contig, config, reference, variants, variant_sample)
                producer_backend = cache.backend
                projections = {h: _projection(root, contig, h, cache.length) for h in ('HAP1', 'HAP2')}
                pending = []
                destinations = {}
                for hap in ('HAP1', 'HAP2'):
                    chosen = [r for r in compact_rows if r['haplotype'] == hap]
                    g = group.create_group(hap)
                    offsets = np.arange(len(chosen)+1, dtype=np.int64)*window_size
                    g.create_dataset('offsets', data=offsets, compressor=codec)
                    g.create_dataset('lengths', data=np.full(len(chosen), window_size, dtype=np.int32), compressor=codec)
                    destinations[hap] = {}
                    for name, dtype in ARRAYS.items():
                        shape = (int(offsets[-1]), len(CHANNELS)) if name == 'features' else (int(offsets[-1]),)
                        # 같은 채널을 이어 압축하고 window 단위로 읽어 작은 파일과 재해독을 줄임.
                        block = min(window_size, 1 << 20)
                        chunks = (block, len(CHANNELS)) if name == 'features' else (block,)
                        destinations[hap][name] = g.create_dataset(name, shape=shape, chunks=chunks, dtype=dtype,
                                                                   order='F' if name == 'features' else 'C', compressor=codec,
                                                                   write_empty_chunks=False)
                    for i, row in enumerate(chosen):
                        pending.append((int(row['raw_start']), hap, i, int(row['actual_length'])))
                # 양 HAP의 인접 reference block을 재사용함. 전체 염색체 배열은 만들지 않음.
                for start, hap, index, actual in sorted(pending):
                    values = projected(cache, projections[hap], start, actual, window_size)
                    for name, value in values.items():
                        destinations[hap][name][index*window_size:(index+1)*window_size] = value
                # 독립 검증의 배열·해독기와 생산 버퍼가 동시에 남지 않도록 반환함.
                cache.close()
                projections.clear()
                if pending:
                    del values, value
            if identities != _inputs(root, contig, alignment, reference, alignment_reference, config, variants):
                raise ValueError('pileup 계산 중 입력이 변경됨')
            hashes = {}
            for hap in ('HAP1', 'HAP2'):
                hashes[hap] = {}
                for name in ARRAYS:
                    a = group[hap][name];h = hashlib.sha256()
                    block = a.chunks[0]
                    for start in range(0, a.shape[0], block):
                        h.update(np.ascontiguousarray(a[start:start+block]))
                    hashes[hap][name] = h.hexdigest()
            metadata['pileup'] = dict(fingerprint_sha256=fingerprint, provenance=provenance,
                                     reference=contract, content_sha256=hashes)
            atomic_json(directory/'contig_metadata.json', metadata)
            os.rename(temporary, directory/'pileup.zarr')
        except BaseException:
            if temporary.exists():
                shutil.rmtree(temporary)
            raise
    print(f'[PILEUP] sample={root.name} contig={contig} phase=source_trace seconds={time.monotonic()-started:.3f} producer={producer_backend}', flush=True)
    validation = {}
    validate_contig(directory, source=True, metrics=validation)
    print(f'[PILEUP] sample={root.name} contig={contig} phase=complete seconds={time.monotonic()-started:.3f} validator={validation["source_backend"]}', flush=True)
    atomic_json(directory/'filter_qc.json', dict(wgs_payload=MODE, filter_policy=config.to_dict(),
                window_index_selection='unchanged', feature_policy=POLICY))
    return dict(record, reused=producer_backend=='reused', production_backend=producer_backend, source_validation=validation)


def validate_contig(directory, *, source=False, metrics=None):
    """각 window를 한 번 읽으며 SHA·구조·source 대조를 모두 수행함."""
    started = time.monotonic()
    directory = Path(directory)
    metadata = json.loads((directory/'contig_metadata.json').read_text())
    sequence = sequence_only.validate_contig(directory, source=source, payload=MODE)
    group = zarr.open_group(str(directory/'pileup.zarr'), mode='r')
    feature_set = group.attrs.get('feature_set')
    channels = resolve_channels(group.attrs.get('channels'), feature_set)
    policy = dict(POLICY, feature_set=feature_set)
    insertion_channels = [i for i, name in enumerate(channels) if name.endswith('insertion_ratio')]
    source_order = [CHANNELS.index(name) for name in channels]
    if (group.attrs.get('schema') != SCHEMA or group.attrs.get('sample') != metadata['sample_id']
            or group.attrs.get('contig') != metadata['contig']
            or group.attrs.get('policy') != policy or group.attrs.get('fingerprint_sha256') != metadata['pileup']['fingerprint_sha256']):
        raise ValueError('pileup schema/채널/정책/지문 불일치')
    source_config = metadata['pileup']['provenance']
    if (source_config.get('schema')!=SCHEMA or source_config.get('channels')!=list(channels)
            or source_config.get('producer')!=feature_set or source_config.get('policy')!=policy
            or digest(source_config)!=metadata['pileup']['fingerprint_sha256']):
        raise ValueError('pileup provenance 계약/지문 불일치')
    bam = cache = None
    try:
        if source:
            identities = source_config['inputs']
            if identities != [file_identity(p[0]) for p in identities]:
                raise ValueError('pileup source 지문 변경')
            bam = pysam.AlignmentFile(identities[0][0], 'r', reference_filename=identities[3][0])
            config = ExportConfig(**source_config['config'])
            reference_contract(bam, identities[1][0], identities[3][0], metadata['contig'])
            hap_reference_contract(directory.parents[2], bam, identities[1][0], metadata['contig'])
            cache = ReferenceCache(bam, metadata['contig'], config, identities[1][0], identities[-2][0], source_config['variant_sample'], oracle=True, block=131072)
        seq = zarr.open_group(str(directory/'hap_sequence.zarr'), mode='r')
        pending = []
        projections, views, hashes = {}, {}, {}
        monotonic = True
        for hap in ('HAP1', 'HAP2'):
            a, h = group[hap], seq[hap]
            n = h['lengths'].shape[0]
            if set(a.array_keys()) != set(ARRAYS) | {'offsets', 'lengths'}:
                raise ValueError('pileup 배열 목록 불일치')
            offsets = a['offsets'][:]; lengths = a['lengths'][:]
            if offsets.shape != (n+1,) or not np.array_equal(lengths, h['lengths'][:]) or offsets[0] != 0 or not np.array_equal(np.diff(offsets), lengths):
                raise ValueError('pileup/window offset 불일치')
            if source:
                projections[hap] = _projection(directory.parents[2], metadata['contig'], hap, cache.length)
            hashes[hap] = {}
            for name, dtype in ARRAYS.items():
                arr = a[name]
                shape = (int(offsets[-1]), len(CHANNELS)) if name == 'features' else (int(offsets[-1]),)
                if arr.dtype != dtype or arr.shape != shape:
                    raise ValueError('pileup dtype/shape 불일치: '+name)
                hashes[hap][name] = hashlib.sha256()
            starts = h['raw_start'][:]; actual = h['actual_length'][:]
            monotonic = monotonic and bool(np.all(starts[1:] >= starts[:-1]))
            views[hap] = (a, offsets)
            pending.extend((int(starts[i]), hap, i, int(actual[i]), int(lengths[i])) for i in range(n))
        # SHA는 저장 순서로 갱신함. 보통의 정렬된 창은 양 HAP의 reference cache도 공유함.
        if source and monotonic:
            pending.sort()
        values_read = 0
        for start, hap, i, actual, padded in pending:
            a, offsets = views[hap]
            observed = {name:np.ascontiguousarray(a[name][offsets[i]:offsets[i+1]]) for name in ARRAYS}
            for name, values in observed.items():
                if name == 'features' and not np.isfinite(values).all():
                    raise ValueError('pileup NaN/inf')
                hashes[hap][name].update(values)
                values_read += values.nbytes
            kinds = observed['mapping_kind']
            if np.any((kinds[:actual] != 1) & (kinds[:actual] != 2)) or np.any(kinds[actual:] != 0):
                raise ValueError('pileup 관측/삽입/padding 상태 오류')
            missing = kinds != 1
            if np.any(observed['reference_base_code'][missing] != 255):
                raise ValueError('비대응 위치에 reference 염기 코드가 있음')
            for name in [key for key in ARRAYS if key not in {'features','mapping_kind','reference_base_code'}]:
                if np.any(observed[name][missing] != 0):
                    raise ValueError('삽입/padding에 가짜 per-base 관측이 있음')
            if np.any(observed['features'][missing][:, [k for k in range(len(channels)) if k not in insertion_channels]] != 0) or np.any(observed['features'][actual:] != 0):
                raise ValueError('비대응 feature 정책 위반')
            if source:
                expected = projected(cache, projections[hap], start, actual, padded)
                if channels != CHANNELS:
                    expected['features'] = expected['features'][:, source_order]
                for name in ARRAYS:
                    if not np.array_equal(observed[name], expected[name]):
                        raise ValueError('독립 reference source 대조 실패: '+name)
        for hap in ('HAP1', 'HAP2'):
            for name, state in hashes[hap].items():
                if state.hexdigest() != metadata['pileup']['content_sha256'][hap][name]:
                    raise ValueError('pileup 내용 SHA 불일치: '+name)
        if source and identities != [file_identity(p[0]) for p in identities]:
            raise ValueError('독립 대조 중 pileup source 지문 변경')
        if metrics is not None:
            metrics.update(seconds=time.monotonic()-started, windows=len(pending),
                           pileup_array_bytes_read=values_read, source_backend=cache.backend if cache else 'not_requested',
                           reference_blocks=cache.blocks_computed if cache else 0, cache_hits=cache.cache_hits if cache else 0)
        return sequence
    finally:
        if cache is not None:
            cache.close()
        if bam is not None:
            bam.close()




def export(root, alignment, reference, alignment_reference, *, contigs, read_type, variants, variant_sample=None, min_mapq=30, jobs=None):
    root = Path(root).resolve()
    variants = Path(variants).resolve()
    variant_sample = variant_sample or root.name
    probe = VariantLookup(variants, variant_sample)
    probe.close()
    sequence_only.preflight(root)
    rows = read_windows(root/'windows/windows.tsv')
    settings = json.loads((root/'windows/windows.json').read_text())
    if set(contigs) != {r['contig'] for r in rows}:
        raise ValueError('전체 contig/window 계약 불일치')
    layout = root/'compact'/LAYOUT
    if any(p.name != LAYOUT for p in (root/'compact').glob('*') if p.is_dir()):
        raise ValueError('다른 payload와 같은 출력 root를 사용하면 안 됨')
    if (layout/'manifest.json').is_file():
        previous = json.loads((layout/'manifest.json').read_text())
        if (previous.get('schema') != SCHEMA or
                resolve_channels(previous.get('channels'), previous.get('feature_set')) != CHANNELS):
            raise ValueError('기존 pileup 채널 순서와 새 출력 계약이 다름; 별도 output/state가 필요함')
    layout.mkdir(parents=True, exist_ok=True)
    manifest = dict(schema=SCHEMA, sample=root.name, wgs_payload=MODE, status='running',
                    validation_level='full', channels=CHANNELS, feature_policy=POLICY, feature_set=FEATURE_SET, records={})
    atomic_json(layout/'manifest.json', manifest)
    cpu = min(int(os.environ.get('MFM_SAMPLE_CPU_BUDGET', '1')), int(os.environ.get('MFM_CPU_POOL_GRANTED', os.environ.get('MFM_SAMPLE_CPU_BUDGET', '1'))))
    workers = min(len(contigs), max(1, cpu-1), jobs or max(1, cpu-1))
    work = [(str(root), str(layout), c, [r for r in rows if r['contig'] == c and r['kept'].lower() in {'1', 'true'}],
             int(settings['window_bp']), str(alignment), str(reference), str(alignment_reference), read_type, min_mapq, str(variants), variant_sample) for c in contigs]
    lease = os.environ.get('MFM_CPU_POOL_LEASE_FILE')
    memory_root = Path(lease).parent.parent if lease else None
    worker_memory, address_space = _worker_memory_limits(memory_active(memory_root))
    for _, record, error in cpu_results(_contig, work, workers, max_workers=workers, parent_cpus=int(cpu > 1),
                                       memory_per_worker_bytes=worker_memory, address_space_per_worker_bytes=address_space,
                                       mp_context=multiprocessing.get_context('spawn')):
        if error:
            raise error
        manifest['records'][record['contig']] = record
        atomic_json(layout/'manifest.json', manifest)
    annotation_records = [dict(contig=c, **{h+'_chain': str(root/'individual_fasta'/h.upper()/(c+'.chain')) for h in ('hap1', 'hap2')}) for c in contigs]
    persist_window_annotations(root, layout, root/'windows', annotation_records, dict(path=str(reference)))
    validate_window_annotations(root)
    manifest['status'] = 'complete'
    atomic_json(layout/'manifest.json', manifest)
    atomic_json(root/'qc/stage12_trace'/MODE/'trace_manifest.json', dict(schema='mfm_stage70_trace_manifest_v1',
                sample=root.name, compact_root=str(layout), wgs_payload=MODE, status='pass',
                contigs={c: 'pass' for c in contigs}, scope='FASTA/window 전수 대조 및 독립 HTSlib reference 집계 전수 대조',
                source_validation={c:manifest['records'][c]['source_validation'] for c in contigs}, completed_at=time.time()))
    stamp(root, layout)
    return manifest


def stamp(root, layout):
    meta_path = root/'metadata/sample.json'
    if os.environ.get('PIPELINE_MANIFEST'):
        import csv
        with Path(os.environ['PIPELINE_MANIFEST']).open() as f:
            matches = [r for r in csv.DictReader(f, delimiter='\t') if r['sample_id'] == root.name]
        if len(matches) != 1:
            raise ValueError('pileup sample manifest 불일치')
        meta = matches[0]
        shutil.copyfile(os.environ['PIPELINE_MANIFEST'], root/'metadata/source_manifest.tsv')
    else:
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    meta.update(sample_id=root.name, sample=root.name, has_rna=bool(meta.get('rna_path') and meta.get('rna_path') != '-'), wgs_payload=MODE, wgs_readset_available=False,
                wgs_depth_available=True, wgs_depth_coordinate_system='reference_projected_to_hap',
                pileup_feature_available=True, pileup_feature_set=FEATURE_SET, pileup_feature_schema=SCHEMA, pileup_feature_channels=list(CHANNELS), compact_root=str(layout), rna_payload=os.environ.get('MFM_RNA_OUTPUT_MODE', 'full'))
    if os.environ.get('MFM_SEQUENCE_DEFER_RNA') == '1' and meta['has_rna']:
        meta.update(completion_scope='dna_pileup_feature_only_rna_deferred', rna_processing_status='deferred',
                    rna_coverage_available=False, ase_available=False, cnv_available=False)
    atomic_json(meta_path, meta)
    atomic_json(root/'metadata/preprocess.json', dict(schema=SCHEMA, sample=root.name, wgs_payload=MODE,
                compact_root=str(layout), channels=CHANNELS, feature_set=FEATURE_SET, feature_policy=POLICY, compressor='zstd', compression_level=5))


