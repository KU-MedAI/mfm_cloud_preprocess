"""한 번 읽기 native 집계기의 고정 count ABI를 연결함."""
import ctypes
import os
import numpy as np
from ..record_stats_parallel import NativeUnavailable, native_tool
from .schema import COUNT_FIELDS, SUM_FIELDS


class Output(ctypes.Structure):
    _fields_ = [('counts', ctypes.c_void_p), ('sums', ctypes.c_void_p)]


class CompleteNative:
    def __init__(self, bam, contig, config, *, oracle=False):
        self.handle = None
        tool = 'pileup_oracle_complete' if oracle else 'pileup_counts_complete'
        if not config.exclude_unmapped:
            raise NativeUnavailable('complete native는 mapped-only 계약을 사용함')
        try:
            self.library = ctypes.CDLL(str(native_tool(tool)))
        except OSError as error:
            raise NativeUnavailable(str(error)) from error
        error_args = [ctypes.c_char_p, ctypes.c_size_t]
        self.library.mfm_complete_open.argtypes = [ctypes.c_char_p]*3+[ctypes.c_int]*3+error_args
        self.library.mfm_complete_open.restype = ctypes.c_void_p
        self.library.mfm_complete_count.argtypes = [ctypes.c_void_p, ctypes.c_int64, ctypes.c_int64, ctypes.POINTER(Output)]+error_args
        self.library.mfm_complete_count.restype = ctypes.c_int
        self.library.mfm_complete_close.argtypes = [ctypes.c_void_p]
        self.library.mfm_complete_close.restype = None
        flags = sum(value for name, value in (('unmapped',4),('secondary',256),('qcfail',512),('duplicate',1024),('supplementary',2048))
                    if getattr(config, 'exclude_'+name))
        error = ctypes.create_string_buffer(2048)
        reference = getattr(bam, 'reference_filename', None)
        self.handle = self.library.mfm_complete_open(os.fsencode(bam.filename), os.fsencode(reference) if reference else b'',
                    contig.encode(), config.min_mapping_quality, flags, int(config.require_proper_pair), error, len(error))
        if not self.handle:
            raise RuntimeError('complete native 초기화 실패: '+error.value.decode(errors='replace'))

    def calculate(self, start, end):
        if not self.handle:
            raise RuntimeError('종료된 complete native')
        counts = np.zeros((3, len(COUNT_FIELDS), end-start), dtype=np.uint32)
        sums = np.zeros((3, len(SUM_FIELDS), end-start), dtype=np.uint64)
        output = Output(counts.ctypes.data, sums.ctypes.data)
        error = ctypes.create_string_buffer(2048)
        if self.library.mfm_complete_count(self.handle, start, end, ctypes.byref(output), error, len(error)):
            raise RuntimeError('complete native 집계 실패: '+error.value.decode(errors='replace'))
        return counts, sums

    def close(self):
        if self.handle:
            self.library.mfm_complete_close(self.handle)
            self.handle = None

    def __del__(self):
        self.close()
