"""성공 증명과 파일 지문이 일치하는 name-collated 준비 BAM만 재사용함."""
from dataclasses import replace
from pathlib import Path
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import time

from .config import PreparedPairedInput

SCHEMA = 'mfm_prepared_bam_cache_v1'


def fingerprint(path):
    path = Path(path)
    stat = path.stat()
    return {'path': str(path.resolve()), 'device': stat.st_dev, 'inode': stat.st_ino,
            'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns, 'ctime_ns': stat.st_ctime_ns}


def implementation_digest():
    base = Path(__file__).resolve().parents[1]
    sources = [base / 'personal_realign.py', *sorted((base / 'realign').glob('*.py'))]
    digest = hashlib.sha256()
    for path in sources:
        digest.update(str(path.relative_to(base)).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


IMPLEMENTATION = implementation_digest()


def cache_path(output):
    return Path(str(output) + '.prepared.json')


def atomic_json(path, value):
    with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False, encoding='utf-8') as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
        temporary = Path(handle.name)
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def key(config, sample, hap, alignment, contigs, route_label):
    source = Path(alignment)
    indexes = [Path(str(source) + suffix) for suffix in ('.bai', '.crai')]
    indexes += [source.with_suffix('.bai'), source.with_suffix('.crai')]
    reference = []
    if source.suffix.lower() == '.cram' and config.reference_fasta:
        for path in [Path(config.reference_fasta), Path(str(config.reference_fasta) + '.fai')]:
            if path.exists():
                reference.append(fingerprint(path))
    executable = shutil.which(config.samtools_bin)
    tool = fingerprint(executable) if executable else {'requested': config.samtools_bin}
    return {'schema': SCHEMA, 'sample': sample,
            'hap': None if config.read_routing_policy == 'all_to_both' else hap,
            'route_label': route_label, 'source': fingerprint(source),
            'indexes': [fingerprint(path) for path in dict.fromkeys(indexes) if path.is_file()],
            'reference': reference, 'contigs': list(contigs), 'samtools': tool,
            'implementation': IMPLEMENTATION,
            'policy': {name: getattr(config, name) for name in (
                'read_routing_policy', 'unknown_read_policy', 'include_unmapped',
                'illumina_duplicate_mode', 'verify_counts', 'count_tolerance_abs')},
            'hap_read_mode': config.hap_read_mode}


def complete_metrics(counts, stats):
    required_counts = {'hp1', 'hp2', 'untagged', 'all_primary_mapped', 'unmapped'}
    required_stats = {'routed_primary', 'rescue_unmapped_primary', 'paired_fragments',
                      'alignment_input_primary', 'alignment_input_routed_primary',
                      'alignment_input_rescue_primary', 'alignment_input_duplicate'}
    if not isinstance(counts, dict) or not isinstance(stats, dict):
        return False
    if not required_counts <= set(counts) or not required_stats <= set(stats):
        return False
    if any(type(counts[k]) is not int or counts[k] < 0 for k in required_counts):
        return False
    if any(type(stats[k]) is not int or stats[k] < 0 for k in required_stats):
        return False
    primary = stats['alignment_input_primary']
    return (primary == 2 * stats['paired_fragments']
            and primary == stats['alignment_input_routed_primary'] + stats['alignment_input_rescue_primary']
            and stats['alignment_input_duplicate'] <= primary)


def load(config, output, expected, input_counts):
    proof = cache_path(output)
    try:
        record = json.loads(proof.read_text())
        if record.get('schema') != SCHEMA or record.get('status') != 'complete':
            return None
        if record.get('key') != expected or record.get('output') != fingerprint(output):
            return None
        counts = record['input_counts']
        stats = record['stats']
        if not complete_metrics(counts, stats):
            return None
        if input_counts is not None and counts != input_counts:
            return None
        if record.get('name_collated') is not True:
            return None
        checked = subprocess.run([config.samtools_bin, 'quickcheck', '-v', str(output)],
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if checked.returncode or record['output'] != fingerprint(output):
            return None
        stats = dict(stats)
        previous = stats.get('input_metrics_execution')
        if previous:
            stats['input_metrics_execution'] = {'mode': 'cached', 'cache_reused': True,
                                               'original_execution': previous}
        stats['prepared_cache'] = {'status': 'reused', 'proof': str(proof),
                                   'created_at': record['created_at']}
        return PreparedPairedInput(Path(output), record['route_label'], stats, counts, True)
    except (OSError, ValueError, KeyError, TypeError):
        return None


def prepare(config, sample, hap, alignment, contigs, work_dir, input_counts,
            route_label, output, builder):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    proof = cache_path(output)
    requested = os.environ.get('MFM_PREPARED_BAM_CACHE', '1')
    if requested not in {'0', '1'}:
        raise ValueError('MFM_PREPARED_BAM_CACHE는0 또는1이어야 함')
    enabled = requested == '1'
    with Path(str(proof) + '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        before = key(config, sample, hap, alignment, contigs, route_label)
        if enabled and not config.overwrite:
            reused = load(config, output, before, input_counts)
            if reused is not None and before == key(config, sample, hap, alignment, contigs, route_label):
                print(f'[PREPARED-CACHE] reuse sample={sample} route={route_label}', flush=True)
                return reused
        proof.unlink(missing_ok=True)
        result = builder(config, sample, hap, alignment, contigs, work_dir, input_counts, route_label)
        if before != key(config, sample, hap, alignment, contigs, route_label):
            raise RuntimeError('준비 BAM 생성 도중 원본/정책 지문이 변경됨. 성공 증명을 저장하지 않음.')
        if enabled:
            if not result.name_collated or not complete_metrics(result.input_counts, result.stats):
                raise RuntimeError('준비 BAM의 성공 계약이 불완전함')
            checked = subprocess.run([config.samtools_bin, 'quickcheck', '-v', str(output)],
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
            if checked.returncode:
                raise RuntimeError('준비 BAM quickcheck 실패. 캐시 성공 증명을 저장하지 않음.')
            created = time.time()
            atomic_json(proof, {'schema': SCHEMA, 'status': 'complete', 'created_at': created,
                               'key': before, 'output': fingerprint(output), 'name_collated': True,
                               'route_label': result.route_label, 'stats': result.stats,
                               'input_counts': result.input_counts})
            stats = dict(result.stats)
            stats['prepared_cache'] = {'status': 'created', 'proof': str(proof), 'created_at': created}
            result = replace(result, stats=stats)
            print(f'[PREPARED-CACHE] created sample={sample} route={route_label}', flush=True)
        return result


def discard(output):
    cache_path(output).unlink(missing_ok=True)
