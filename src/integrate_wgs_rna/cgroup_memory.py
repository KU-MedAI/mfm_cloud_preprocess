"""비활성 파일 캐시만 보수적으로 회수 가능하다고 보는 입장 RAM 판정."""
from pathlib import Path


def pressure_usage(current, maximum, folder=None):
    if current is None or not maximum:
        return current
    try:
        if folder is None:
            relative = next(line.split('::', 1)[1].strip().lstrip('/')
                            for line in Path('/proc/self/cgroup').read_text().splitlines() if '::' in line)
            mount = Path('/sys/fs/cgroup')
            leaf = mount / relative
            folder = next(parent for parent in [leaf, *leaf.parents]
                          if (parent == mount or mount in parent.parents)
                          and (parent / 'memory.max').is_file()
                          and (parent / 'memory.max').read_text().strip() == str(maximum))
        actual = Path(folder) / 'memory.current'
        if actual.is_file():
            current = max(int(current), int(actual.read_text()))
        values = dict(line.split() for line in (Path(folder) / 'memory.stat').read_text().splitlines())
        # shmem·dirty·writeback 및 활성 캐시를 제외하고 비활성 파일의10%도 여유로 남긴다.
        clean_file = max(0, int(values['file']) - int(values.get('shmem', 0))
                         - int(values.get('file_dirty', 0)) - int(values.get('file_writeback', 0)))
        clean_inactive = max(0, int(values['inactive_file']) - int(values.get('file_dirty', 0))
                             - int(values.get('file_writeback', 0)))
        reclaimable = min(int(current), clean_file, clean_inactive)
        return int(current) - reclaimable * 9 // 10
    except (OSError, ValueError, KeyError, StopIteration):
        return current
