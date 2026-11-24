"""HTSlib 검증기에 bounded count 배열만 전달함. Python 생산기와 독립된 집계 경로임."""
import ctypes
import os
from ..record_stats_parallel import NativeUnavailable, native_tool

FIELDS = ('bases','mapq','baseq','baseq_n','forward','deletion','insertion',
          'deletion_event','max_deletion','supplementary')


class Output(ctypes.Structure):
    _fields_ = [(name,ctypes.c_void_p) for name in FIELDS]


class NativeOracle:
    tool = 'pileup_oracle'

    def __init__(self, bam, contig, config):
        self.handle = None
        if not config.exclude_unmapped:
            raise NativeUnavailable('이 native pileup은 mapped-only 정책에서 사용함')
        try:
            self.library = ctypes.CDLL(str(native_tool(self.tool)))
        except OSError as error:
            raise NativeUnavailable('native pileup 라이브러리 로드 불가: '+str(error)) from error
        error_args = [ctypes.c_char_p,ctypes.c_size_t]
        self.library.mfm_pileup_open.argtypes = [ctypes.c_char_p]*3+[ctypes.c_int]*3+error_args
        self.library.mfm_pileup_open.restype = ctypes.c_void_p
        self.library.mfm_pileup_count.argtypes = [ctypes.c_void_p,ctypes.c_int64,ctypes.c_int64,ctypes.POINTER(Output)]+error_args
        self.library.mfm_pileup_count.restype = ctypes.c_int
        self.library.mfm_pileup_close.argtypes = [ctypes.c_void_p]
        self.library.mfm_pileup_close.restype = None
        flags = sum(flag for name,flag in (('unmapped',4),('secondary',256),('qcfail',512),('duplicate',1024),('supplementary',2048))
                    if getattr(config,'exclude_'+name))
        error = ctypes.create_string_buffer(2048)
        reference = getattr(bam,'reference_filename',None)
        self.handle = self.library.mfm_pileup_open(os.fsencode(bam.filename), os.fsencode(reference) if reference else b'',
                contig.encode(), config.min_mapping_quality, flags, int(config.require_proper_pair), error, len(error))
        if not self.handle:
            raise RuntimeError(self.tool+' 입력 초기화 실패: '+error.value.decode(errors='replace'))

    def calculate(self, start, end):
        from .compute import Counts
        if self.handle is None:
            raise RuntimeError('종료된 native oracle')
        counts = Counts(start,end)
        output = Output(*(getattr(counts,name).ctypes.data for name in FIELDS))
        error = ctypes.create_string_buffer(2048)
        if self.library.mfm_pileup_count(self.handle,start,end,ctypes.byref(output),error,len(error)):
            raise RuntimeError(self.tool+' 계산 실패: '+error.value.decode(errors='replace'))
        return counts.finish()

    def close(self):
        if self.handle:
            self.library.mfm_pileup_close(self.handle)
            self.handle = None

    def __del__(self):
        self.close()


class NativeProducer(NativeOracle):
    tool = 'pileup_counts'
