#!/usr/bin/env python3
"""파일 inventory와 해시를 확인하고 실행 데이터 없이 클라우드 코드를 내보냄."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = {
    'docker/Dockerfile', '.dockerignore', 'docker/config.pileup.example.json',
    'docker/entrypoint.sh', 'docker/smoke.py', 'docs/Docker_실행.md',
    'run_pipeline.sh', 'config_v2.sh', 'handoff_config.example.json',
    'cloud_config.pileup.example.json', 'environment.yml', 'requirements.txt',
    'SOURCE_MANIFEST.json', 'GUIDE.md', 'docs/클라우드_실행.md',
    'manifests/공통자산.json',
    'references/hg38_autosomes_blacklist_union.bed',
    'src/integrate_wgs_rna/pileup_features/export.py',
    'src/integrate_wgs_rna/native/pileup_counts.cpp',
    'src/integrate_wgs_rna/native/pileup_oracle.cpp',
}


def package_file(root: Path, relative: str) -> Path:
    """inventory의 경로 이탈과 외부/내부 symlink를 모두 거부함."""
    rel = Path(relative)
    if rel.is_absolute() or not rel.parts or '..' in rel.parts:
        raise ValueError('허용되지 않은 파일 경로: ' + relative)
    path = root / rel
    if root.resolve() not in path.resolve().parents:
        raise ValueError('패키지 밖 파일 경로: ' + relative)
    for parent in (path, *path.parents):
        if parent == root:
            break
        if parent.is_symlink():
            raise ValueError('패키지 파일은 symlink일 수 없음: ' + relative)
    return path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 << 20), b''):
            value.update(block)
    return value.hexdigest()


def verify_package(root: Path) -> tuple[dict, list[dict]]:
    root = root.resolve()
    errors, records = [], []
    manifest = root / 'CLOUD_MANIFEST.json'
    try:
        data = json.loads(manifest.read_text())
        if not isinstance(data, dict) or data.get('schema') != 'mfm_cloud_package_v1':
            raise ValueError('클라우드 inventory schema 불일치')
        pending = data['files']
        if not isinstance(pending, list) or any(
            not isinstance(row, dict) or not isinstance(row.get('path'), str)
            or type(row.get('size_bytes')) is not int or row['size_bytes'] < 0
            or not isinstance(row.get('sha256'), str) for row in pending
        ):
            raise ValueError('클라우드 inventory 파일 레코드 형식 오류')
        records = pending
        names = [row['path'] for row in records]
        if len(names) != len(set(names)):
            raise ValueError('클라우드 inventory 경로 중복')
        for name in sorted(REQUIRED - set(names)):
            errors.append('필수 파일 inventory 누락: ' + name)
        for row in records:
            path = package_file(root, row['path'])
            if not path.is_file():
                errors.append('파일 누락: ' + row['path'])
            elif path.stat().st_size != row['size_bytes'] or digest(path) != row['sha256']:
                errors.append('파일 크기/해시 불일치: ' + row['path'])
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append(str(error))
    return {
        'status': 'pass' if not errors else 'fail',
        'files_checked': len(records),
        'code_and_metadata_bytes': sum(row.get('size_bytes', 0) for row in records),
        'errors': errors,
        'scope': '코드·문서·입력목록·작은 blacklist 주석; 원본/완료데이터·클라우드 실행 검증 아님',
    }, records


def export_package(root: Path, destination: Path) -> dict:
    report, records = verify_package(root)
    if report['errors']:
        raise ValueError('검증 실패 패키지는 내보내지 않음: ' + '; '.join(report['errors'][:3]))
    destination = destination.expanduser().absolute()
    if destination.exists() or destination.is_symlink():
        raise ValueError('내보낼 새 경로가 이미 존재함: ' + str(destination))
    if root.resolve() == destination.resolve() or root.resolve() in destination.resolve().parents:
        raise ValueError('패키지 내부로 내보낼 수 없음')
    destination.mkdir(parents=True)
    for row in records:
        target = destination / row['path']
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(package_file(root, row['path']), target)
    shutil.copy2(root / 'CLOUD_MANIFEST.json', destination / 'CLOUD_MANIFEST.json')
    copied, _ = verify_package(destination)
    if copied['errors']:
        raise ValueError('복사본 검증 실패: ' + '; '.join(copied['errors'][:3]))
    return dict(copied, exported_to=str(destination), runtime_data_exported=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--export-to', type=Path, help='존재하지 않는 새 경로로 코드만 내보냄')
    args = parser.parse_args()
    try:
        report = export_package(ROOT, args.export_to) if args.export_to else verify_package(ROOT)[0]
    except (OSError, ValueError) as error:
        report = {'status': 'fail', 'errors': [str(error)]}
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return int(report['status'] != 'pass')


if __name__ == '__main__':
    raise SystemExit(main())
