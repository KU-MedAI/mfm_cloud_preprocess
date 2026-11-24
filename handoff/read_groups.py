"""분담/배치를 보존한 리드별 목록을 기존 handoff 실행으로 연결함."""
from pathlib import Path
from types import SimpleNamespace
import argparse
import copy
import fcntl
import json
import subprocess
import sys
import tempfile

from . import cli, planner
from .queue_resources import effective_memory, full_ram

ROOT = Path(__file__).resolve().parents[1]


def plans(lab, read_type, batch, config):
    directory = ROOT / 'manifests/리드별' / lab
    batches = sorted(p.name for p in directory.iterdir() if p.is_dir()) if batch == 'all' else [batch.zfill(3)]
    results = []
    for number in batches:
        selection = directory / number / (read_type + '.tsv')
        if not selection.is_file():
            raise ValueError('해당 기관/배치별 목록이 없음: ' + str(selection))
        rows = cli.selection(SimpleNamespace(selection=str(selection)))
        expected_lab = 'TGIL' if lab == 'tgil' else 'MedAI'
        if any(r['assigned_lab'] != expected_lab or r['batch_id'] != number or r['read_type'] != read_type for r in rows):
            raise ValueError('리드 분리 목록과 기존 분담 배정이 다름')
        limit = effective_memory(config) * 9 // 10
        blocked = {r['sample_id']: '공통16GiB·다운로드12GB 포함 보수 RAM예약 초과: '
                   + str(full_ram([r], config) + 12_000_000_000) + '/' + str(limit)
                   for r in rows if full_ram([r], config) + 12_000_000_000 > limit}
        eligible = [r for r in rows if r['sample_id'] not in blocked]
        # 순서는 분담 TSV의 원래 순서를 유지하고 Scheduler에서 우선순위를 적용함.
        item_config = copy.deepcopy(config)
        item_config.pop('assignment_snapshot', None)  # 신규 논리배치는 현행 master를 검증함.
        item_config['state_root'] = str(Path(config['state_root']) / 'read_groups' / read_type)
        item_config.setdefault('pipeline_env', {}).setdefault('MFM_CPU_POOL_DIR',
            str(Path(config['state_root']) / 'cpu_pool' / 'read_groups'))
        item_config['queue_policy'] = {'mode': 'long_read_first' if read_type == 'long_read' else 'ready_first'}
        results.append(dict(lab=lab, batch=number, read_type=read_type,
                            selection=str(selection), rows=eligible, blocked=blocked,
                            total=len(rows), config=item_config))
    return results


def summary(item):
    return {k: item[k] for k in ['lab', 'batch', 'read_type', 'total', 'blocked']} | {'eligible': len(item['rows'])}


def dry_run(items):
    results = []
    # 임시 경로에 config/selection만 작성하고 데이터/output/state 경로는 만들지 않음.
    with tempfile.TemporaryDirectory(prefix='mfm_read_group_') as temporary:
        for item in items:
            result = summary(item)
            if item['rows']:
                config = Path(temporary) / 'config.json'
                selected = Path(temporary) / 'selected.tsv'
                config.write_text(json.dumps(item['config']))
                planner.write_tsv(selected, item['rows'], fields=list(item['rows'][0]))
                process = subprocess.run([sys.executable, '-B', '-m', 'handoff.cli', 'dry-run',
                    '--selection', str(selected), '--config', str(config), '--executor', 'local'],
                    cwd=ROOT, capture_output=True, text=True)
                if process.returncode:
                    raise RuntimeError('리드별 기존 진입점 dry-run 실패: ' + process.stderr[-2000:])
                result['pipeline_plan'] = json.loads(process.stdout)
            results.append(result)
    print(json.dumps({'status': 'dry_run', 'batches': results,
        'selected': sum(x['eligible'] for x in results),
        'blocked': sum(len(x['blocked']) for x in results),
        'side_effects': '데이터 다운로드·설치·전처리·state 생성 없음. 임시 설정/선택 파일만 검사 후 정리함.'},
        ensure_ascii=False, indent=2))
    return 0


def run(items, config):
    state = Path(config['state_root']) / 'read_groups'
    state.mkdir(parents=True, exist_ok=True)
    with (state / 'controller.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        results = []
        for item in items:
            directory = state / item['read_type'] / (item['lab'] + '_' + item['batch']) / 'launch'
            directory.mkdir(parents=True, exist_ok=True)
            result = summary(item)
            selected = directory / 'selected.tsv'
            if selected.exists() and planner.read_tsv(selected) != item['rows']:
                raise ValueError('기존 실행 목록과 RAM 통과 선택이 달라짐. 별도 state_root로 실행해야 함.')
            from .local_queue import locked, read
            old_config = read(directory / 'config.json')
            if old_config and old_config != item['config'] and any(
                    locked(directory.parent, r['sample_id']) for r in item['rows']):
                raise ValueError('실행 중 worker의 config를 덮어쓸 수 없음')
            cli.atomic_json(directory / '제외사유.json', item['blocked'])
            cli.atomic_json(directory / 'config.json', item['config'])
            if item['rows']:
                selected = directory / 'selected.tsv'
                temporary = directory / '.selected.tsv.tmp'
                planner.write_tsv(temporary, item['rows'], fields=list(item['rows'][0]))
                temporary.replace(selected)
                command = ['bash', str(ROOT / 'run_pipeline.sh'), 'run', '--selection', str(selected),
                           '--config', str(directory / 'config.json'), '--executor', 'local']
                result['returncode'] = subprocess.run(command, cwd=ROOT).returncode
            else:
                result['returncode'] = 0
            queue_status = read(directory.parent / 'queue/status.json')
            if result['returncode'] == 0 and queue_status.get('target_reached'):
                result['target_reached'] = True
                results.append(result)
                cli.atomic_json(state / (item['lab'] + '_' + item['read_type'] + '_summary.json'),
                    {'batches': results, 'status': 'target_reached', 'blocked': sum(len(x['blocked']) for x in results)})
                return 0
            results.append(result)
            cli.atomic_json(state / (item['lab'] + '_' + item['read_type'] + '_summary.json'),
                {'batches': results, 'blocked': sum(len(x['blocked']) for x in results),
                 'status': 'failed' if result['returncode'] else 'running'})
            if result['returncode']:
                return result['returncode']
        pending = sum(len(x['blocked']) for x in results)
        report = {'batches': results, 'blocked': pending,
                  'status': 'complete_with_pending' if pending else 'complete'}
        cli.atomic_json(state / (items[0]['lab'] + '_' + items[0]['read_type'] + '_summary.json'), report)
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 1 if pending else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['dry-run', 'run'])
    parser.add_argument('--lab', required=True, choices=['medai', 'tgil'])
    parser.add_argument('--read-type', required=True, choices=['long_read', 'short_read'])
    parser.add_argument('--batch', default='all', help='기존 배치 번호 또는 all')
    parser.add_argument('--config', default=str(ROOT / 'handoff_config.example.json'))
    args = parser.parse_args()
    config_path = Path(args.config)
    if not config_path.is_file():
        raise ValueError('설정 파일이 없음: ' + str(config_path))
    config = cli.load_config(config_path)
    items = plans(args.lab, args.read_type, args.batch, config)
    return dry_run(items) if args.command == 'dry-run' else run(items, config)


if __name__ == '__main__':
    raise SystemExit(main())
