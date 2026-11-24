"""RNA 전체 해독 검사와 손상 입력 보존·대체 다운로드를 담당한다."""
from __future__ import annotations
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import time
import uuid
from pathlib import Path
from .downloads import fetch, fetch_or_index, safe_path, verify
from .input_check import fingerprint
from .resume import sample_downloads

SCHEMA = 'mfm_rna_full_decode_v1'

class InputCorruptionError(ValueError):
    """내용 손상이 확인된 입력에만 새 다운로드를 허용한다."""

def atomic_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False, encoding='utf-8') as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
        temporary = stream.name
    os.replace(temporary, path)

def validate_full(path, report, sample, samtools='samtools', threads=4, decode_reference=None):
    """헤더/EOF와 전체 레코드 해독을 검사하고 동일 파일·검사기 결과만 재사용한다."""
    path = Path(path)
    report = Path(report)
    if path.suffix == '.cram':
        from .io_integrity import validate
        from .alignment_integrity import AlignmentCorruption
        try:
            validate(path,Path(str(report)+'.physical.json'))
        except AlignmentCorruption as error:
            raise InputCorruptionError('RNA CRAM 전체 물리 블록 손상: '+str(error)) from error
    before = fingerprint(path)
    executable = shutil.which(str(samtools))
    if not executable:
        raise FileNotFoundError('samtools 실행 파일 없음: ' + str(samtools))
    validator = fingerprint(executable)
    reference = fingerprint(decode_reference) if decode_reference and path.suffix == '.cram' else None
    identity = {'schema':SCHEMA, 'sample':sample, 'input':before,
                'validator':validator, 'decode_reference':reference}
    if report.is_file() and time.time_ns() - max(before['mtime_ns'], before['ctime_ns']) > 2_000_000_000:
        try:
            previous = json.loads(report.read_text())
            if previous.get('status') == 'pass' and all(previous.get(key) == value for key,value in identity.items()):
                return previous
        except (OSError, ValueError, TypeError):
            pass
    payload = dict(identity, status='checking', started=time.time())
    commands = [[executable,'quickcheck','-v',str(path)],
                [executable,'flagstat','-@',str(max(1,int(threads))),'-O','json']]
    if reference:
        commands[1] += ['--input-fmt-option', 'reference=' + reference['path']]
    commands[1].append(str(path))
    for label,command in zip(['quickcheck','full_decode'],commands):
        result = subprocess.run(command,capture_output=True,text=True)
        if result.returncode:
            detail = (result.stdout + result.stderr)[-8000:]
            payload.update(status='failed',failed_check=label,returncode=result.returncode,
                           error=detail,finished=time.time())
            atomic_json(report,payload)
            if before != fingerprint(path):
                raise RuntimeError('RNA 검사 중 입력이 변경되었습니다: ' + str(path))
            # signal 종료/OOM/명령 실패는 입력 손상으로 분류하지 않는다.
            markers = ('bgzf','truncated','error reading','could not read header',
                       'not in sequence format','not sequence data','eof marker is absent')
            system_errors = ('permission denied','cannot allocate memory','out of memory',
                             'memory allocation','too many open files','no space left',
                             'input/output error','failed to populate reference','failed to fetch reference')
            if result.returncode > 0 and not any(mark in detail.lower() for mark in system_errors) and (label == 'quickcheck' or any(mark in detail.lower() for mark in markers)):
                raise InputCorruptionError(f'RNA {label} 입력 손상: {path}: {detail}')
            raise RuntimeError(f'RNA {label} 검사 명령 실패: {path}: {detail}')
        if label == 'full_decode':
            try:
                flagstat = json.loads(result.stdout)
                if not isinstance(flagstat,dict) or not isinstance(flagstat.get('QC-passed reads'),dict):
                    raise ValueError('flagstat JSON 계약 불일치')
            except ValueError as error:
                raise RuntimeError('RNA 전체 해독 결과 JSON 오류: ' + str(path)) from error
            payload['flagstat'] = flagstat
    if before != fingerprint(path):
        raise RuntimeError('RNA 검사 중 입력이 변경되었습니다: ' + str(path))
    payload.update(status='pass',finished=time.time())
    atomic_json(report,payload)
    return payload

def prepare_rna(config, row, resources, decode_reference=None):
    """같은 URL의 RNA만 최대2회 별도 경로에 확보하고 검증한 경로를 worker에 전달한다."""
    from .sequence_mode import rna_required
    if not rna_required(config,row):
        return config, resources
    alignments = [resource for resource in resources if resource['kind'] == 'rna_alignment']
    if len(alignments) != 1:
        raise ValueError('pair 샘플 RNA alignment resource는 정확히1개여야 합니다: ' + row['sample_id'])
    original = dict(alignments[0])
    sample = row['sample_id']
    data = Path(config['data_root'])
    initial = safe_path(data,original['path'])
    if not initial.is_file():
        fetch(original,data,config['state_root'])
    st = initial.stat()
    source_fingerprint = {'path':str(initial.resolve()),'size':st.st_size,
                          'mtime_ns':st.st_mtime_ns,'ctime_ns':st.st_ctime_ns,'inode':st.st_ino}
    key = hashlib.sha256((sample+'\n'+original['path']+'\n'+original['url']).encode()).hexdigest()
    mapping = data/'.rna_recovery_records'/(key+'.json')
    recovery_events = data/'.rna_recovery_records'/(key+'.events.jsonl')
    local = dict(config)
    current = resources
    if mapping.is_file():
        try:
            previous = json.loads(mapping.read_text())
            if (previous.get('schema') == 'mfm_rna_recovery_binding_v1'
                and previous.get('source_fingerprint') == source_fingerprint
                and previous.get('sample') == sample and previous.get('url') == original['url']):
                relative = previous['replacement_relative']
                candidate = safe_path(data,relative)
                if candidate.is_file() and verify(candidate,original.get('size_bytes'),original.get('checksum'),original.get('checksum_type','md5')):
                    local['sample_rna_paths'] = dict(config.get('sample_rna_paths',{}), **{sample:relative})
                    current = sample_downloads(local,row,resources)
        except (OSError, ValueError, KeyError, TypeError):
            pass
    attempts = config.get('rna_recovery_attempts',2)
    if isinstance(attempts,bool) or not isinstance(attempts,int) or not 0 <= attempts <= 5:
        raise ValueError('rna_recovery_attempts는0~5 정수여야 합니다.')
    for attempt in range(attempts + 1):
        alignment = next(resource for resource in current if resource['kind'] == 'rna_alignment')
        path = safe_path(data,alignment['path'])
        report_key = hashlib.sha256(str(path).encode()).hexdigest()
        report = Path(config['state_root'])/'rna_input_validation'/sample/(report_key+'.json')
        try:
            if not verify(path,alignment.get('size_bytes'),alignment.get('checksum'),alignment.get('checksum_type','md5')):
                atomic_json(report,{'schema':SCHEMA,'sample':sample,'status':'failed',
                    'failed_check':'size_or_checksum','path':str(path),'finished':time.time()})
                raise InputCorruptionError('RNA 크기/checksum 검증 실패: ' + str(path))
            validate_full(path,report,sample,config.get('samtools','samtools'),
                          min(4,config['sample_cpus']),decode_reference)
            for resource in current:
                if resource['kind'] == 'rna_index':
                    fetch_or_index(resource,data,config['state_root'],min(4,config['sample_cpus']))
            if alignment['path'] != original['path']:
                atomic_json(mapping,{'schema':'mfm_rna_recovery_binding_v1','sample':sample,
                    'url':original['url'],'source_fingerprint':source_fingerprint,
                    'replacement_relative':alignment['path'],'validation_report':str(report),
                    'validated_at':time.time()})
            local['rna_full_decode_report'] = str(report)
            return local,current
        except InputCorruptionError as error:
            recovery_events.parent.mkdir(parents=True,exist_ok=True)
            with recovery_events.open('a') as stream:
                stream.write(json.dumps({'sample':sample,'time':time.time(),'failed_path':str(path),
                    'validation_report':str(report),'attempt':attempt,'error':str(error)},ensure_ascii=False)+'\n')
            if attempt == attempts:
                raise InputCorruptionError(f'RNA 손상 복구 {attempts}회 후 실패. 입력/QC 기록 보존: {sample}') from error
            relative = str(Path('rna_recovery')/sample/uuid.uuid4().hex/Path(original['path']).name)
            local['sample_rna_paths'] = dict(local.get('sample_rna_paths',{}), **{sample:relative})
            current = sample_downloads(local,row,resources)
            replacement = next(resource for resource in current if resource['kind'] == 'rna_alignment')
            print(f'[RNA-RECOVERY] sample={sample} retry={attempt+1}/{attempts} preserved={path} replacement={relative}',flush=True)
            fetch(replacement,data,config['state_root'])
