"""명시적 복구에서만 성공한 HiFi 정렬 산출물을 확인하고 후처리부터 이어감."""
import hashlib
import json
import os
from pathlib import Path
import shlex


def completed_output_proof(config, sample, hap, alignment, target, output, log_dir):
    flag = os.environ.get('MFM_REALIGN_RESUME_COMPLETED_OUTPUT', '0')
    if flag not in {'0', '1'}:
        raise ValueError('MFM_REALIGN_RESUME_COMPLETED_OUTPUT은0/1이어야 함')
    if flag == '0' or not output.exists():
        return None
    if (config.overwrite or config.dry_run or str(config.platform).lower() not in {'hifi', 'pacbio_hifi'}
            or os.environ.get('MFM_WGS_BACKEND', 'cpu') != 'cpu'):
        raise ValueError('완료 산출물 복구는 overwrite/dry-run 없는 CPU HiFi에서만 허용함')
    script = log_dir / f'{sample}.HAP{hap}.realign.sh'
    resource = log_dir / f'{sample}.HAP{hap}.resource.json'
    usage = json.loads(resource.read_text())
    if usage.get('exit_code') != 0 or usage.get('command') != ['bash', str(script)]:
        raise ValueError('완료 정렬의 성공 명령 증거가 일치하지 않음')
    assignments = {}
    text = script.read_text()
    for line in text.splitlines():
        if '=' in line and line.split('=', 1)[0] in {'INPUT', 'TARGET', 'OUTPUT', 'PART'}:
            key, value = line.split('=', 1)
            parts = shlex.split(value)
            if len(parts) != 1:
                raise ValueError('정렬 script 경로가 단일 값이 아님')
            assignments[key] = parts[0]
    expected = dict(INPUT=alignment, TARGET=target, OUTPUT=output, PART=Path(str(output)+'.part'))
    if any(key not in assignments or Path(assignments[key]).resolve() != path.resolve()
           for key, path in expected.items()):
        raise ValueError('정렬 script의 원본/참조/출력 경로가 현재 실행과 다름')
    if '-ax map-hifi' not in text or 'mv "$PART" "$OUTPUT"' not in text:
        raise ValueError('검증된 HiFi 원자 완료 script가 아님')
    if expected['PART'].exists() or output.stat().st_size < 1:
        raise ValueError('불완전한 정렬 산출물은 복구하지 않음')
    started = float(usage['started_at_unix'])
    finished = started + float(usage['elapsed_sec'])
    if any(path.stat().st_mtime > started+1 for path in [alignment, target, script]):
        raise ValueError('정렬 시작 이후 원본/참조/script가 변경됨')
    if not started-1 <= output.stat().st_mtime <= finished+2:
        raise ValueError('완료 CRAM의 생성 시각이 성공 명령 범위와 다름')
    return dict(status='verified_completed_alignment', script=str(script),
                script_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),
                resource_usage=str(resource), output=str(output), output_bytes=output.stat().st_size,
                output_mtime=output.stat().st_mtime, original_alignment_seconds=usage['elapsed_sec'],
                original_cpu_seconds=usage['user_cpu_sec']+usage['system_cpu_sec'],
                next_actions='기존 native index·전수 read/HP/PS/QC·개인축 검증을 실행함')
