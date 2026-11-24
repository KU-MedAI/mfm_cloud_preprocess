"""운영중선택은SHA로고정한이전배정에대조하고신규선택은현행배정에대조함."""
from pathlib import Path
import hashlib,json

def master_path(root,args):
    default=Path(root)/'manifests/분담/전체_분담배치.tsv'
    config=getattr(args,'effective_config',None)
    if not config:return default
    spec=json.loads(Path(config).read_text()).get('assignment_snapshot')
    if spec is None:return default
    if not isinstance(spec,dict) or not isinstance(spec.get('relative_path'),str):raise ValueError('운영 배정 snapshot 형식 오류')
    p=(Path(root)/spec['relative_path']).resolve();directory=(Path(root)/'manifests/분담').resolve()
    if directory not in p.parents or p.is_symlink():raise ValueError('운영 배정 snapshot은 패키지 분담 경로여야 함')
    if hashlib.sha256(p.read_bytes()).hexdigest()!=spec.get('sha256'):raise ValueError('운영 배정 snapshot SHA 불일치')
    return p
