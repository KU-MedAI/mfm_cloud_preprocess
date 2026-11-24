"""동봉 manifest를 이용해 공통 자산을 다운로드하거나 원본에서 생성한다."""
from __future__ import annotations

import argparse
import fcntl
import gzip
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

from .bundled_assets import get as bundled_get, manifest_path
from .downloads import fetch, safe_path, verify

NCBI = 'references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
BROAD = 'references/Homo_sapiens_assembly38.fasta'
AUTOSOMES = 'references/hg38_autosomes.fa'
HGSVC_AUTOSOMES = 'references/hgsvc_hg38_autosomes.fa'
RNA_EXTRA = 'references/rna_index_nonautosome_gtex.fa'


def dependencies(relative):
    if relative.endswith('.fai'):
        return [relative[:-4]]
    if relative == AUTOSOMES:
        return [NCBI + '.fai']
    if relative == RNA_EXTRA:
        return [BROAD + '.fai']
    return []


def operation(relative):
    if relative.endswith('.fai'):
        return 'samtools faidx로 생성'
    if relative == NCBI:
        return '원본 gzip 다운로드·해제'
    if relative == AUTOSOMES:
        return 'NCBI chr1–22 추출; FASTA 행 너비 60'
    if relative == HGSVC_AUTOSOMES:
        return 'HGSVC 공식 no-ALT 원본 MD5 검증·chr1–22 대문자 추출; 행 너비 50'
    if relative == RNA_EXTRA:
        return 'Broad 비상염색체 추출; alt/decoy/HLA 제외; 행 너비 100'
    return '공개 URL 다운로드'


def ordered(records, selected):
    result = []
    visited = set()
    visiting = set()

    def add(relative):
        if relative in visited:
            return
        if relative in visiting or relative not in records:
            raise ValueError('자산 dependency 오류: ' + relative)
        visiting.add(relative)
        for name in dependencies(relative):
            add(name)
        visiting.remove(relative)
        visited.add(relative)
        result.append(relative)

    for name in selected:
        add(name)
    return result


def subset(source, target, names, width, uppercase=False):
    """FASTA를 고정 순서·행 너비로 쓴다. 염색체 전체를 메모리에 올리지 않는다."""
    import pysam

    block_size = width * 10000
    with pysam.FastaFile(str(source)) as fasta, target.open('w') as handle:
        for name in names:
            handle.write('>' + name + '\n')
            length = fasta.get_reference_length(name)
            for start in range(0, length, block_size):
                sequence = fasta.fetch(name, start, min(length, start + block_size))
                if uppercase:
                    sequence = sequence.upper()
                handle.write('\n'.join(sequence[i:i + width] for i in range(0, len(sequence), width)) + '\n')


class Restorer:
    def __init__(self, root, base, manifest, samtools):
        self.root = root
        self.base = base
        self.manifest = manifest
        self.records = {r['path']: r for r in manifest['files']}
        self.samtools = samtools
        # 기존 bundled 검증 cache를 사용한다. standalone 준비용이며 cleanup을 호출하지 않는다.
        self.config = {'bundled_assets_root': str(base), 'data_root': str(base)}
        self.done = {}
        self.history = []

    def finish(self, relative, temporary):
        row = self.records[relative]
        if not verify(temporary, row['size_bytes'], row['sha256'], 'sha256'):
            raise ValueError('생성 파일 size/SHA256 불일치; .part 보존: ' + relative)
        path = safe_path(self.base, relative)
        os.replace(temporary, path)
        return path

    def restore(self, relative):
        if relative in self.done:
            return self.done[relative]
        row = self.records[relative]
        path = safe_path(self.base, relative)
        for name in dependencies(relative):
            self.restore(name)
        if path.exists() or path.is_symlink():
            # 기존 파일이 다르면 명확히 중단하며 덮어쓰지 않는다.
            path = bundled_get(self.config, self.root, relative, required=True)
            self.history.append({'path': relative, 'status': 'reused'})
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_name(path.name + '.part')
            if temporary.is_symlink():
                raise ValueError('임시 파일이 symlink임: ' + str(temporary))
            if relative.endswith('.fai'):
                subprocess.run([self.samtools, 'faidx', '--fai-idx', str(temporary),
                                str(self.done[relative[:-4]])], check=True)
                path = self.finish(relative, temporary)
            elif relative == AUTOSOMES:
                subset(self.done[NCBI], temporary, [f'chr{i}' for i in range(1, 23)], 60)
                path = self.finish(relative, temporary)
            elif relative == HGSVC_AUTOSOMES:
                from .references import HGSVC, decompress, index
                compressed = fetch({'url':HGSVC,'path':'.download_cache/hgsvc_hg38.no_alt.fa.gz',
                                    'size_bytes':963047116,'checksum_type':'md5',
                                    'checksum':'6a5ceceed91ae5c14026454f09c4b6c4'}, self.base, self.base)
                source = self.base/'.download_cache/hgsvc_hg38.no_alt.fa'
                source.parent.mkdir(parents=True,exist_ok=True)
                source = decompress(compressed, source)
                index(source, self.samtools)
                subset(source, temporary, [f'chr{i}' for i in range(1,23)], 50, uppercase=True)
                path = self.finish(relative, temporary)
            elif relative == RNA_EXTRA:
                import pysam
                from rna_profile import select_nonautosome
                with pysam.FastaFile(str(self.done[BROAD])) as fasta:
                    names = select_nonautosome(list(fasta.references))
                subset(self.done[BROAD], temporary, names, 100)
                path = self.finish(relative, temporary)
            elif relative == NCBI:
                compressed = fetch({'url': row['source_url'],
                                    'path': '.download_cache/' + path.name + '.gz',
                                    'size_bytes': 872949833, 'checksum_type': 'md5',
                                    'checksum': self.manifest['source_validation']['ncbi_gzip_md5']},
                                   self.base, self.base)
                with gzip.open(compressed, 'rb') as source, temporary.open('wb') as output:
                    shutil.copyfileobj(source, output, 8 << 20)
                path = self.finish(relative, temporary)
            else:
                if not row.get('source_url'):
                    raise ValueError('공개 URL/생성 규칙 없음: ' + relative)
                path = fetch({'url': row['source_url'], 'path': relative,
                              'size_bytes': row['size_bytes'], 'checksum_type': 'sha256',
                              'checksum': row['sha256']}, self.base, self.base)
            self.history.append({'path': relative, 'status': 'restored', 'operation': operation(relative)})
        self.done[relative] = path
        print(f"[자산 준비] {relative}: {self.history[-1]['status']}", flush=True)
        return path


def main(argv=None):
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assets-root', type=Path, default=root / 'assets',
                        help='manifest.json이 있는 자산 폴더; 기본값 ./assets')
    parser.add_argument('--group', choices=['all', 'references', 'panel'], default='all')
    parser.add_argument('--file', action='append', default=[],
                        help='manifest 상대 경로를 지정; 여러 번 가능하며 --group 대신 적용')
    parser.add_argument('--dry-run', action='store_true', help='계획만 출력; 파일/폴더 생성·다운로드 없음')
    parser.add_argument('--samtools', default=os.environ.get('SAMTOOLS', 'samtools'))
    args = parser.parse_args(argv)
    base = args.assets_root.resolve()
    manifest = json.loads(manifest_path({'bundled_assets_root': str(base)}, root).read_text())
    if manifest.get('schema') != 'mfm_bundled_assets_v1':
        raise ValueError('assets manifest schema 불일치')
    records = {r['path']: r for r in manifest['files']}
    if len(records) != len(manifest['files']):
        raise ValueError('manifest 경로 중복')
    for relative, row in records.items():
        safe_path(base, relative)
        if not relative.startswith(('references/', 'common_vcf/')):
            raise ValueError('지원하지 않는 자산 경로: ' + relative)
        if int(row['size_bytes']) < 1 or len(row['sha256']) != 64:
            raise ValueError('자산 크기/SHA256 누락: ' + relative)
    selected = args.file or [name for name in records if args.group == 'all'
                            or name.startswith('references/' if args.group == 'references' else 'common_vcf/')]
    sequence = ordered(records, selected)
    plan = {'operation': '공통 자산 준비', 'assets_root': str(base),
            'selected_files': len(selected), 'files_with_dependencies': len(sequence),
            'final_selected_bytes': sum(records[name]['size_bytes'] for name in sequence),
            'files': [{'path': name, 'operation': operation(name),
                       'source_url': records[name].get('source_url'), 'dependencies': dependencies(name),
                       'size_bytes': records[name]['size_bytes'], 'sha256': records[name]['sha256'],
                       'present': (base / name).is_file()} for name in sequence]}
    if args.dry_run:
        print(json.dumps(plan, ensure_ascii=False, indent=2))
        return 0
    if any(name.endswith('.fai') for name in sequence) and shutil.which(args.samtools) is None:
        raise ValueError('samtools 필요; jj_mfm을 활성화하거나 --samtools 지정')
    if AUTOSOMES in sequence or RNA_EXTRA in sequence:
        try:
            import pysam
        except ImportError as error:
            raise ValueError('파생 FASTA 생성에는 pysam 필요; jj_mfm 환경 사용') from error
    if shutil.which('curl') is None:
        raise ValueError('curl 필요')
    restorer = Restorer(root, base, manifest, args.samtools)
    report = {'schema': 'mfm_asset_restore_v1', 'status': 'running', 'started': time.time(),
              'group': args.group, 'files': restorer.history}
    with (base / '.asset_download.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            for relative in sequence:
                restorer.restore(relative)
            report['status'] = 'pass'
        except Exception as error:
            report.update(status='fail', error=str(error))
            raise
        finally:
            report['finished'] = time.time()
            temporary = base / '.asset_restore_report.json.part'
            temporary.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
            os.replace(temporary, base / '.asset_restore_report.json')
    print(f"[완료] {len(restorer.done)}개 자산 size/SHA256 준비 확인")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
