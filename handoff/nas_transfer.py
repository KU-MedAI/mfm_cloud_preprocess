"""완료 샘플의 수동 NAS 배치 준비·검증 확인·대용량 로컬 자료 회수를 수행함."""
from contextlib import contextmanager
from pathlib import Path, PurePosixPath
import argparse
import csv
import fcntl
import importlib.util
import json
import os
import re
import secrets
import shutil
import sys
import time

from . import archive_receipt

ROOT = Path(__file__).resolve().parents[1]
PAYLOADS = {'compact', 'rna_personal'}


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    with temporary.open('w') as stream:
        stream.write(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


@contextmanager
def lock(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a+') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        yield


def safe_relative(name):
    relative = PurePosixPath(name)
    if (not name or relative.is_absolute() or '..' in relative.parts or '\n' in name
            or '\r' in name or '\x00' in name or str(relative) != name):
        raise ValueError('허용되지 않은 상대 경로')
    return relative


def files(root):
    root = Path(root)
    def walk(directory):
        for path in sorted(directory.iterdir()):
            relative = path.relative_to(root).as_posix()
            if relative == archive_receipt.RECEIPT:
                continue
            if path.is_symlink():
                raise ValueError('이관 패키지 symlink 거부: ' + str(path))
            if path.is_dir():
                yield from walk(path)
            elif path.is_file():
                yield path, relative
            else:
                raise ValueError('일반 파일이 아닌 이관 항목: ' + str(path))
    yield from walk(root)


def catalog_snapshot(root, project):
    name = '27_make_catalog.py' if (ROOT / '09_personal_fasta.sh').is_file() else '25_make_catalog.py'
    spec = importlib.util.spec_from_file_location('nas_catalog_snapshot', ROOT / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.scan_sample(root, project, {}, {})


def inventories(root, output):
    count = total = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + '.tmp')
    with temporary.open('w') as stream:
        for path, relative in files(root):
            before = path.stat()
            sha = archive_receipt.digest(path)
            after = path.stat()
            if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
                    after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
                raise ValueError('목록 생성 중 파일 변경: ' + relative)
            stream.write(json.dumps(dict(path=relative, size=after.st_size, sha256=sha),
                                    ensure_ascii=False) + '\n')
            count += 1
            total += after.st_size
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, output)
    return dict(inventory_path=str(output), inventory_sha256=archive_receipt.digest(output),
                file_count=count, total_bytes=total)


def batch_path(config, batch_id):
    if not re.fullmatch(r'[A-Za-z0-9_-]{1,80}', batch_id):
        raise ValueError('잘못된 이관 배치 ID')
    return Path(config['state_root']) / 'nas_exports' / batch_id / 'batch.json'


def checked_local(config, row, root):
    from . import cli
    from .sequence_mode import completion_mode, payload_matches
    from .resume import sample_environment
    if (not cli.processing_complete(root,row['sample_id'],completion_mode(config,row),require_cleanup=True)
            or not payload_matches(config,row,root)):
        return False
    expected = sample_environment(config,row['sample_id']).get('WINDOW_STORAGE_POLICY')
    policy = json.loads((root/'metadata/window_policy.json').read_text())
    return expected is None or policy.get('storage_policy') == expected


def prepare(config, rows, destination_root, batch_size=10, allow_partial=False, batch_id=None):
    from . import cli
    from .sequence_mode import completion_mode
    from .local_queue import locked
    if not 1 <= batch_size <= 10:
        raise ValueError('수동 이관 배치는 1~10명이어야 함')
    destination_root = Path(destination_root)
    if (not destination_root.is_absolute() or destination_root.name != 'final_pileup'
            or '..' in destination_root.parts):
        raise ValueError('NAS의 절대 final_pileup 경로를 지정해야 함')
    directory = Path(config['state_root']) / 'nas_exports'
    def same_context(value):
        if (value.get('schema')!='mfm_nas_batch_v1' or value.get('run_id')!=config['run_id']
                or value.get('selection_fingerprint')!=config['selection_fingerprint']
                or value.get('destination_root')!=str(destination_root)):
            raise ValueError('기존 이관 배치의 실행/선택/NAS 경로와 다름')
        return value
    with lock(directory / 'prepare.lock'):
        if batch_id:
            value = json.loads(batch_path(config, batch_id).read_text())
            return same_context(value)
        for path in sorted(directory.glob('*/batch.json')):
            value = json.loads(path.read_text())
            if value['status'] == 'prepared':
                return same_context(value)
        ready = []
        for row in rows:
            sample = row['sample_id']
            root = cli.sample_root(config, row)
            record = Path(config['state_root']) / 'samples' / (sample + '.json')
            if not record.is_file() or json.loads(record.read_text()).get('status') != 'success':
                continue
            if archive_receipt.matches(config, row, root) is not None or locked(config['state_root'], sample):
                continue
            if not checked_local(config,row,root):
                continue
            ready.append(row)
        if not ready or (len(ready) < batch_size and not allow_partial):
            return dict(status='waiting',ready=len(ready),required=batch_size)
        selected = ready[:batch_size]
        identity = time.strftime('%Y%m%d_%H%M%S') + '_' + secrets.token_hex(4)
        path = batch_path(config, identity)
        value = dict(schema='mfm_nas_batch_v1',status='preparing',batch_id=identity,
                     token=secrets.token_hex(24),run_id=config['run_id'],
                     selection_fingerprint=config['selection_fingerprint'],
                     destination_root=str(destination_root),created=time.time(),samples=[])
        for row in selected:
            sample = row['sample_id']
            if not re.fullmatch(r'[A-Za-z0-9_-]+', sample):
                raise ValueError('허용되지 않은 sample ID')
            root = cli.sample_root(config,row)
            if root.is_symlink() or not (root / 'compact').is_dir():
                raise ValueError('완료 payload가 없음')
            root = root.resolve()
            if any(path.name not in {'compact','rna_personal','qc','metadata','logs'} for path in root.iterdir()):
                raise ValueError('모델 전용 정리 후의 패키지만 이관 가능함')
            relative = Path(config['output_root']).name + '/' + '/'.join(
                [row['data_group'],row['project_id'],sample])
            safe_relative(relative)
            item = dict(sample=sample,source_root=str(root),destination=str(destination_root / relative),
                        mode=completion_mode(config,row),relative_destination=relative,
                        model_bytes=json.loads((Path(config['state_root'])/'samples'/(sample+'.json')).read_text())['model_bytes'],
                        catalog_snapshot=catalog_snapshot(root,row['project_id']))
            item.update(inventories(root,path.parent/(sample+'.files.jsonl')))
            if locked(config['state_root'],sample):
                raise ValueError('목록 생성 중 샘플이 실행됨')
            item['retained_hashes'] = {name:archive_receipt.digest(root/name)
                                       for name in sorted(archive_receipt.RETAINED)}
            value['samples'].append(item)
        value['status']='prepared'
        save(path,value)
        return value


def verify_source(root, inventory, allow_missing=False):
    expected={}
    with Path(inventory).open() as stream:
        for line in stream:
            row=json.loads(line);safe_relative(row['path']);expected[row['path']]=row
    observed=set()
    for path,relative in files(root):
        row=expected.get(relative)
        if row is None or path.stat().st_size!=row['size'] or archive_receipt.digest(path)!=row['sha256']:
            raise ValueError('NAS 검증 뒤 source 파일 변경: '+relative)
        observed.add(relative)
    missing=set(expected)-observed
    if missing and (not allow_missing or any(PurePosixPath(name).parts[0] not in PAYLOADS for name in missing)):
        raise ValueError('검증한 source 파일 누락')
    return observed


def commit(config, rows, receipt, remove_source=False):
    from . import cli
    from .local_queue import locked
    path=batch_path(config,receipt.get('batch_id',''))
    with lock(path.parent/'commit.lock'):
        batch=json.loads(path.read_text())
        if (receipt.get('token')!=batch['token'] or receipt.get('verification')!='sha256_all_files'
                or batch['run_id']!=config['run_id']
                or batch['selection_fingerprint']!=config['selection_fingerprint']):
            raise ValueError('다른 실행/배치의 NAS 검증 결과임')
        item=next((r for r in batch['samples'] if r['sample']==receipt.get('sample')),None)
        if item is None or any(receipt.get(k)!=item[k] for k in [
                'destination','inventory_sha256','file_count','total_bytes']):
            raise ValueError('NAS 수신 경로/내용 검증 계약 불일치')
        row=next(r for r in rows if r['sample_id']==item['sample'])
        root=cli.sample_root(config,row).resolve()
        raw_root=cli.sample_root(config,row)
        if (raw_root.is_symlink() or str(root)!=item['source_root']
                or locked(config['state_root'],item['sample'])):
            raise ValueError('현재 실행 중이거나 다른 source 경로임')
        inventory=Path(item['inventory_path'])
        if inventory.parent!=path.parent or archive_receipt.digest(inventory)!=item['inventory_sha256']:
            raise ValueError('서버 파일 목록 변경')
        previous=archive_receipt.read(root)
        if previous is not None and previous['batch_id']!=batch['batch_id']:
            raise ValueError('다른 배치로 이미 이관된 샘플임')
        if previous is None and not checked_local(config,row,root):
            raise ValueError('현재 완료/QC/출력 정책과 준비한 패키지가 다름')
        if previous is not None:
            archive_receipt.matches(config,row,root)
        observed=verify_source(root,inventory,allow_missing=previous is not None)
        record=dict(item,schema=archive_receipt.SCHEMA,status='nas_verified',batch_id=batch['batch_id'],batch_file=str(path),
                    run_id=batch['run_id'],selection_fingerprint=batch['selection_fingerprint'],
                    verification='sha256_all_files',verified_at=receipt.get('verified_at',time.time()),
                    receiver=receipt.get('receiver',''),local_payload_removed=False)
        if previous is not None and previous['status']=='archived':
            record=previous
        save(root/archive_receipt.RECEIPT,record)
        if remove_source:
            from .completion_cleanup import _opened
            opened=_opened({str(root/name) for name in observed if PurePosixPath(name).parts[0] in PAYLOADS})
            if opened:
                raise ValueError('서버 payload를 읽는 프로세스가 있음: '+str(opened[:3]))
            # payload만 회수하며 이미 실행 중인 구 큐의 QC/완료 계약도 보존함.
            sample_lock=Path(config['state_root'])/'samples'/(item['sample']+'.lock')
            with sample_lock.open('a+') as handle:
                fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
                for name in sorted(PAYLOADS):
                    target=root/name
                    if target.is_symlink():
                        raise ValueError('삭제 payload symlink 거부')
                    if target.exists():
                        shutil.rmtree(target)
                record.update(status='archived',local_payload_removed=True,removed_at=time.time())
                save(root/archive_receipt.RECEIPT,record)
        info_path=Path(config['state_root'])/'samples'/(item['sample']+'.json')
        info=json.loads(info_path.read_text())
        info.update(archive_status=record['status'],artifact_location='nas',artifact_root=item['destination'])
        save(info_path,info)
        item['archive_status']=record['status']
        if all(r.get('archive_status') in {'nas_verified','archived'} for r in batch['samples']):
            batch['status']='complete'
        save(path,batch)
        return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','commit'])
    parser.add_argument('--effective-config',required=True)
    parser.add_argument('--destination-root')
    parser.add_argument('--batch-size',type=int,default=10)
    parser.add_argument('--allow-partial',action='store_true')
    parser.add_argument('--batch')
    parser.add_argument('--remove-source',action='store_true')
    args=parser.parse_args()
    config=json.loads(Path(args.effective_config).read_text())
    with (Path(config['state_root'])/'selected_samples.tsv').open() as stream:
        rows=list(csv.DictReader(stream,delimiter='\t'))
    if args.action=='prepare':
        if not args.destination_root:parser.error('--destination-root 필요함')
        result=prepare(config,rows,args.destination_root,args.batch_size,args.allow_partial,args.batch)
    else:
        result=commit(config,rows,json.load(sys.stdin),args.remove_source)
    print(json.dumps(result,ensure_ascii=False))


if __name__=='__main__':
    main()
