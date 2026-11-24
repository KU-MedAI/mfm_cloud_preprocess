"""실제 실행 cgroup의 CPU·RAM 할당을 읽음."""
from pathlib import Path


def cgroup_limits():
    current = None
    memory = []
    cpus = []
    try:
        relative = next(line.split('::', 1)[1].strip().lstrip('/')
                        for line in Path('/proc/self/cgroup').read_text().splitlines() if '::' in line)
        mount = Path('/sys/fs/cgroup')
        leaf = mount / relative
        if (leaf / 'memory.current').exists():
            current = int((leaf / 'memory.current').read_text())
        for parent in [leaf, *leaf.parents]:
            if parent != mount and mount not in parent.parents:
                break
            value = (parent / 'memory.max').read_text().strip() if (parent / 'memory.max').exists() else 'max'
            if value != 'max':
                memory.append(int(value))
            quota = (parent / 'cpu.max').read_text().split() if (parent / 'cpu.max').exists() else ['max']
            if quota[0] != 'max':
                cpus.append(max(1, int(quota[0]) // int(quota[1])))
    except (OSError, ValueError, StopIteration):
        pass
    return current, min(memory) if memory else None, min(cpus) if cpus else None

