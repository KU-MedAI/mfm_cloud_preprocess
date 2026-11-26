#!/usr/bin/env python3
"""NAS 측에서 완료본을 수동 또는 주기 조회로 pull·검증·선택적 서버 회수함."""
from pathlib import Path, PurePosixPath
import argparse
import fcntl
import hashlib
import json
import os
import shlex
import shutil
import socket
import subprocess
import sys
import time


def log(message):
    print(time.strftime('[%Y-%m-%d %H:%M:%S] ')+message,flush=True)


def digest(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8<<20),b''):value.update(block)
    return value.hexdigest()


def verify(directory, inventory, expected):
    if digest(inventory)!=expected['inventory_sha256']:
        raise ValueError('수신한 서버 파일 목록 SHA 불일치')
    names=set();count=total=0
    with Path(inventory).open() as stream:
        for line in stream:
            row=json.loads(line);name=row['path'];relative=PurePosixPath(name)
            if (relative.is_absolute() or '..' in relative.parts or not relative.parts
                    or str(relative)!=name or name in names):raise ValueError('잘못된 수신 상대 경로')
            path=Path(directory)/relative
            if any(p.is_symlink() for p in [path,*path.parents] if p==Path(directory) or Path(directory) in p.parents):
                raise ValueError('수신 패키지 symlink 거부')
            if not path.is_file() or path.stat().st_size!=row['size'] or digest(path)!=row['sha256']:
                raise ValueError('NAS 파일 checksum 불일치: '+name)
            names.add(name);count+=1;total+=row['size']
    actual=set()
    for base,dirs,files in os.walk(directory,followlinks=False):
        if any((Path(base)/d).is_symlink() for d in dirs):raise ValueError('NAS 디렉터리 symlink 거부')
        for name in files:actual.add((Path(base)/name).relative_to(directory).as_posix())
    if actual!=names or count!=expected['file_count'] or total!=expected['total_bytes']:
        raise ValueError('NAS 파일 개수/경로/용량 불일치')


def remote(args, action, **options):
    command=['env','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1',
             str(Path(args.source_python)),str(Path(args.source_code)/'tools/nas_transfer.py'),
             action,'--effective-config',args.source_config]
    for key,value in options.items():
        if value is True:command.append('--'+key.replace('_','-'))
        elif value is not None and value is not False:command.extend(['--'+key.replace('_','-'),str(value)])
    return [*args.ssh,args.source,shlex.join(command)]


def receive(args,batch):
    root=Path(args.destination).resolve()
    records=root/'_transfer_records'/batch['batch_id'];records.mkdir(parents=True,exist_ok=True)
    progress='--info=stats2' if getattr(args,'watch',False) else '--info=progress2'
    rsync=['rsync','-rltHs','--fsync','--partial-dir=.rsync-partial',progress,'-e',shlex.join(args.ssh)]
    for item in batch['samples']:
        relative=PurePosixPath(item['relative_destination'])
        if relative.is_absolute() or '..' in relative.parts:raise ValueError('수신 경로 이탈')
        final=root/relative
        if str(final)!=item['destination']:raise ValueError('서버와 NAS 수신 위치가 다름')
        inventory=records/(item['sample']+'.files.jsonl')
        subprocess.run([*rsync,'--checksum',args.source+':'+item['inventory_path'],str(inventory)],check=True)
        staging=root/'.incoming'/batch['batch_id']/item['sample']
        for candidate in [final,staging]:
            if any(p.is_symlink() for p in [candidate,*candidate.parents] if p==root or root in p.parents):
                raise ValueError('NAS 수신 경로 symlink 거부')
        if not final.exists():
            staging.mkdir(parents=True,exist_ok=True)
            if shutil.disk_usage(root).free<item['total_bytes']:
                raise ValueError('NAS의 샘플 수신 여유 공간이 부족함')
            log('[NAS] '+item['sample']+' 복사')
            resume_options=['--checksum'] if any(staging.iterdir()) else []
            subprocess.run([*rsync,*resume_options,args.source+':'+item['source_root']+'/',str(staging)+'/'],check=True)
            verify(staging,inventory,item)
            final.parent.mkdir(parents=True,exist_ok=True)
            os.rename(staging,final)
            descriptor=os.open(final.parent,os.O_RDONLY|os.O_DIRECTORY)
            try:os.fsync(descriptor)
            finally:os.close(descriptor)
        else:
            verify(final,inventory,item)
        receipt={key:item[key] for key in ['sample','destination','inventory_sha256','file_count','total_bytes']}
        receipt.update(batch_id=batch['batch_id'],token=batch['token'],verification='sha256_all_files',
                       verified_at=time.time(),receiver=socket.gethostname())
        subprocess.run(remote(args,'commit',remove_source=args.remove_source),input=json.dumps(receipt),
                       text=True,capture_output=True,check=True)
        path=records/(item['sample']+'.verified.json')
        with path.open('w') as stream:
            stream.write(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n');stream.flush();os.fsync(stream.fileno())
        log('[NAS] '+item['sample']+' 검증·이관 완료'+(' / 서버 payload 정리' if args.remove_source else ''))


def pull_once(args):
    response=subprocess.run(remote(args,'prepare',destination_root=str(Path(args.destination).resolve()),
        batch_size=args.batch_size,allow_partial=args.allow_partial or args.watch,batch=args.batch),
        capture_output=True,text=True,check=True)
    batch=json.loads(response.stdout)
    if batch['status']=='waiting':
        log('완료 '+str(batch['ready'])+'명 / '+str(batch['required'])+'명 대기');return False
    if batch['status'] not in {'prepared','complete'}:raise ValueError('수신 가능한 배치 상태가 아님')
    log('이관 배치: '+batch['batch_id'])
    receive(args,batch)
    return True


def run(args):
    while True:
        try:received=pull_once(args)
        except subprocess.CalledProcessError as error:
            if not args.watch:raise
            log('접속/전송 실패. '+str(args.poll_interval)+'초 후 같은 작업을 다시 확인함.')
            if error.stderr:print(error.stderr.strip(),file=sys.stderr,flush=True)
            received=False
        if not args.watch:return
        if not received:time.sleep(args.poll_interval)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True,help='SSH 사용자@원본호스트')
    parser.add_argument('--port',type=int,default=22)
    parser.add_argument('--identity',required=True)
    parser.add_argument('--source-code',required=True)
    parser.add_argument('--source-config',required=True)
    parser.add_argument('--source-python',required=True,help='원본 jj_mfm의 Python 절대 경로')
    parser.add_argument('--destination',required=True)
    parser.add_argument('--batch-size',type=int,default=10)
    parser.add_argument('--allow-partial',action='store_true')
    parser.add_argument('--batch')
    parser.add_argument('--remove-source',action='store_true')
    parser.add_argument('--watch',action='store_true',help='새 완료본을 반복 조회하며 10명 미만도 수신함')
    parser.add_argument('--poll-interval',type=int,default=300,help='대기/접속 실패 후 재조회 간격(초), 기본300')
    args=parser.parse_args()
    if not 1<=args.batch_size<=10:parser.error('--batch-size는1~10이어야 함')
    if args.poll_interval<=0:parser.error('--poll-interval은 양수여야 함')
    if args.watch and args.batch:parser.error('--watch와 특정 --batch는 함께 사용할 수 없음')
    if not shutil.which('ssh') or not shutil.which('rsync'):parser.error('NAS 측 ssh와 rsync 설치가 필요함')
    if args.source.startswith('-') or any(c.isspace() for c in args.source):parser.error('잘못된 SSH 원본')
    root=Path(args.destination)
    if not root.is_absolute() or root.name!='final_pileup' or not root.parent.is_dir():
        parser.error('기존 NAS 상위 경로 아래의 절대 final_pileup 경로가 필요함')
    root.mkdir(exist_ok=True)
    args.ssh=['ssh','-p',str(args.port),'-i',str(Path(args.identity).expanduser()),'-o','IdentitiesOnly=yes',
              '-o','Compression=no']
    if args.watch:
        args.ssh.extend(['-o','BatchMode=yes','-o','ConnectTimeout=15',
                         '-o','ServerAliveInterval=30','-o','ServerAliveCountMax=3'])
    with (root/'.nas_pull.lock').open('a+') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:parser.error('이 NAS 경로에서 이미 수신 도구가 실행 중임')
        if args.watch:log('자동 수신 시작. 새 완료본이 없으면 '+str(args.poll_interval)+'초 후 재조회함.')
        run(args)


if __name__=='__main__':
    try:main()
    except KeyboardInterrupt:
        log('수신 중지. 부분 전송과 이관 기록은 다음 실행에서 재개 가능함.')
        raise SystemExit(130)
    except subprocess.CalledProcessError as error:
        if error.stderr:print(error.stderr, file=sys.stderr)
        raise
