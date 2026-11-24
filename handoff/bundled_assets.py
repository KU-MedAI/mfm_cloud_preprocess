"""동봉한 공통 입력을 검증·재사용하며 raw/cache cleanup과 분리함."""
from __future__ import annotations
import fcntl
import hashlib
import json
import os
import time
from pathlib import Path
from .downloads import safe_path


def asset_root(config, root):
    return Path(config.get('bundled_assets_root', Path(root) / 'assets')).resolve()


def manifest_path(config, root):
    local = asset_root(config, root) / 'manifest.json'
    return local if local.is_file() else Path(root) / 'manifests/공통자산.json'


def inventory(config, root):
    base = asset_root(config, root)
    manifest = manifest_path(config, root)
    if not manifest.is_file():
        return base, {}
    data = json.loads(manifest.read_text())
    if data.get('schema') != 'mfm_bundled_assets_v1':
        raise ValueError('동봉 공통 자산 manifest schema 불일치')
    records = {r['path']: r for r in data['files']}
    if len(records) != len(data['files']):
        raise ValueError('동봉 공통 자산 경로 중복')
    return base, records


def get(config, root, relative, required=False):
    base, records = inventory(config, root)
    record = records.get(relative)
    if record is None:
        if required:
            raise ValueError('동봉 필수 자산 누락: ' + relative)
        return None
    path = safe_path(base, relative)
    if not path.exists() and not path.is_symlink():
        from .asset_restore import Restorer
        base.mkdir(parents=True, exist_ok=True)
        manifest = json.loads(manifest_path(config, root).read_text())
        with (base / '.asset_download.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            Restorer(root, base, manifest, config.get('samtools', 'samtools')).restore(relative)
    if not path.is_file() or path.is_symlink():
        raise ValueError('동봉 자산은 실제 파일이어야 함: ' + str(path))
    cache = Path(config['data_root']) / '.bundled_asset_records'
    cache.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(str(path).encode()).hexdigest()
    report = cache / (key + '.json')
    with (cache / (key + '.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        stat = path.stat()
        signature = {'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
                     'ctime_ns': stat.st_ctime_ns, 'inode': stat.st_ino}
        if stat.st_size != int(record['size_bytes']):
            raise ValueError('동봉 자산 size 불일치: ' + relative)
        previous = json.loads(report.read_text()) if report.is_file() else {}
        # timestamp 해상도가 낮은 파일시스템의 같은 tick 내 덮어쓰기를 보호함.
        recent = time.time_ns() - max(stat.st_mtime_ns, stat.st_ctime_ns) < 2_000_000_000
        if not recent and previous.get('signature') == signature and previous.get('sha256') == record['sha256']:
            return path
        digest = hashlib.sha256()
        with path.open('rb') as handle:
            while block := handle.read(8 << 20):
                digest.update(block)
        if digest.hexdigest() != record['sha256']:
            raise ValueError('동봉 자산 SHA256 불일치: ' + relative)
        after = path.stat()
        if signature != {'size': after.st_size, 'mtime_ns': after.st_mtime_ns,
                         'ctime_ns': after.st_ctime_ns, 'inode': after.st_ino}:
            raise ValueError('동봉 자산 검증 중 파일 변경: ' + relative)
        temporary = report.with_suffix('.part')
        temporary.write_text(json.dumps({'path': str(path), 'sha256': digest.hexdigest(),
                                          'signature': signature}) + '\n')
        os.replace(temporary, report)
    return path


def references(config, root, needs_rna=False, needs_illumina=False, needs_hgsvc=False):
    base, records = inventory(config, root)
    if not records:
        return None
    paths = {
        'ncbi': 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        'autosomes': 'hg38_autosomes.fa',
        'hgsvc': 'hgsvc_hg38_autosomes.fa',
        'broad': 'Homo_sapiens_assembly38.fasta',
        'illumina': 'GRCh38_full_analysis_set_plus_decoy_hla.fa',
        'gtf': 'gencode.v47.annotation.gtf.gz',
        'rna_extra': 'rna_index_nonautosome_gtex.fa',
    }
    needed = ['ncbi', 'autosomes']
    if needs_hgsvc:
        needed += ['hgsvc']
    if needs_rna or needs_illumina:
        needed += ['broad']
    if needs_rna:
        needed += ['gtf', 'rna_extra']
    if needs_illumina:
        needed += ['illumina']
    for role in needed:
        relative = 'references/' + paths[role]
        get(config, root, relative, required=True)
        if role != 'gtf':
            get(config, root, relative + '.fai', required=True)
    return {role: str(base / 'references' / filename) for role, filename in paths.items()}


def panel_directory(config, root):
    base, records = inventory(config, root)
    if any(name.startswith('common_vcf/') for name in records):
        return base / 'common_vcf'
    return Path(config['data_root']) / 'common_vcf'
