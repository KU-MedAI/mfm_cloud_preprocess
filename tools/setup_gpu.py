#!/usr/bin/env python3
"""선택 GPU Python 의존성과 호환 STAR/BWA 도구를 jj_mfm에 준비함."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tarfile,tempfile,urllib.request

ROOT=Path(__file__).resolve().parents[1]
SOURCES={'STAR':('2.7.2a','https://github.com/alexdobin/STAR/archive/refs/tags/2.7.2a.tar.gz'),
         'bwa':('0.7.17','https://github.com/lh3/bwa/archive/refs/tags/v0.7.17.tar.gz')}

def install_tools():
    destination=ROOT/'tools/gpu';destination.mkdir(parents=True,exist_ok=True)
    manifest={}
    for tool,(version,url) in SOURCES.items():
        output=destination/('STAR_2.7.2a' if tool=='STAR' else 'bwa')
        if output.exists():
            if tool=='STAR':
                observed=subprocess.run([str(output),'--version'],text=True,capture_output=True,check=True).stdout.strip()
                if observed!=version:raise ValueError('기존 STAR GPU 호환 도구 버전 불일치')
            manifest[tool]={'version':version,'path':str(output),'status':'reused'};continue
        with tempfile.TemporaryDirectory(prefix='gpu_tool_',dir=destination) as temp:
            directory=Path(temp);archive=directory/'source.tar.gz'
            with urllib.request.urlopen(url,timeout=60) as response,archive.open('wb') as handle:shutil.copyfileobj(response,handle)
            digest=hashlib.sha256(archive.read_bytes()).hexdigest()
            expected={'STAR':'939db81f7b39e2de33815c7760fad1a4e660b6972f641a2e6edcced61a515b8a','bwa':'980b9591b61c60042c4a39b9e31ccaad8d17ff179d44d347997825da3fdf47fd'}
            if digest!=expected[tool]:raise ValueError('GPU 도구 source SHA256 불일치: '+tool)
            with tarfile.open(archive) as source:
                for member in source.getmembers():
                    name=Path(member.name)
                    if name.is_absolute() or '..' in name.parts:raise ValueError('도구 archive 경로 오류')
                source.extractall(directory,filter='data')
            source=next(p for p in directory.iterdir() if p.is_dir())
            if tool=='STAR':binary=source/'bin/Linux_x86_64_static/STAR'
            else:
                subprocess.run(['make','-j','1','CFLAGS=-g -Wall -Wno-unused-function -O2 -fcommon'],cwd=source,check=True)
                binary=source/'bwa'
            if not binary.is_file():raise RuntimeError('GPU 호환 도구 생성 실패: '+tool)
            temporary=output.with_suffix('.part');shutil.copy2(binary,temporary);temporary.chmod(0o755);os.replace(temporary,output)
            manifest[tool]={'version':version,'path':str(output),'source':url,'source_sha256':digest,
                            'binary_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'status':'installed'}
    path=destination/'install_manifest.json'
    old=json.loads(path.read_text()) if path.exists() else {}
    old.update({key:value for key,value in manifest.items() if value['status']!='reused' or key not in old})
    path.write_text(json.dumps(old,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(manifest,ensure_ascii=False,indent=2))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--tools',action='store_true');p.add_argument('--verify',action='store_true')
    a=p.parse_args()
    if a.verify:
        import cupy
        print('CuPy',cupy.__version__,'CUDA runtime',cupy.cuda.runtime.runtimeGetVersion());return
    if Path(sys.prefix).name!='jj_mfm':raise RuntimeError('jj_mfm 환경에서 실행해야 함')
    if a.tools:install_tools();return
    subprocess.run([sys.executable,'-m','pip','install','-r',str(ROOT/'requirements-gpu.txt')],check=True)
    print('GPU 배열 의존성 설치 완료. GPU 정렬은 Parabricks 실행환경과 install-gpu-tools를 추가 준비함.')

if __name__=='__main__':main()
