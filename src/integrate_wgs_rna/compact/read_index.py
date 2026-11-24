"""read 메타정보를 고정 폭으로 저장하고 완료 좌표의 식별 키를 반환한다."""
from array import array


class ReadIndex:
    def __init__(self):
        self.keys = {}
        self.haplotypes = array('B')
        self.phase_sets = array('q')
        self.lengths = array('i')
        self.reference_ends = array('q')
        self.retired = 0
        self.last_boundary = -1

    def retire_before(self, boundary: int) -> None:
        """개인 좌표의 다음 윈도에 닿지 않는 정렬만 반환한다."""
        if boundary < self.last_boundary:
            raise ValueError('개인 윈도 좌표는 증가 순서여야 합니다.')
        self.last_boundary = boundary
        expired = [key for key, index in self.keys.items()
                   if self.reference_ends[index] >= 0 and self.reference_ends[index] < boundary]
        for key in expired:
            del self.keys[key]
        self.retired += len(expired)
