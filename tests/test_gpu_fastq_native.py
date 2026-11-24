"""네이티브 FASTQ의 바이트 보존·오류 검출·FIFO 첫 입력 재생을 대조한다."""
from pathlib import Path
from unittest.mock import patch
import gzip
import io
import os
import tempfile
import unittest
from integrate_wgs_rna.gpu_fastq import split_fastq
from integrate_wgs_rna.gpu_stream_inputs import _Replay
from integrate_wgs_rna.record_stats_parallel import native_tool


class NativeFastqTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        native_tool('gpu_fastq')

    def split(self, raw, paired, mode, replay=False):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ, {'MFM_GPU_FASTQ_BACKEND':mode}):
            paths = [Path(tmp)/f'{i}.gz' for i in range(2 if paired else 1)]
            source = _Replay(raw[:17].splitlines(keepends=True), io.BytesIO(raw[17:])) if replay else io.BytesIO(raw)
            execution = {}
            count = split_fastq(source, paths, paired, cpus=16, producer_cpus=3, execution=execution)
            self.assertEqual(execution['backend'], mode)
            return count, [gzip.decompress(p.read_bytes()) for p in paths]

    def test_태그_CRLF_마지막개행_프레임경계_일치(self):
        for paired in [False, True]:
            records = []
            for i in range(12000):
                for mate in range(2 if paired else 1):
                    records.append(f'@r{i}/{mate+1} HP:i:1\tPS:i:23\tOD:i:1\r\n'.encode()+b'AGCT'*100+b'\r\n+\r\n'+b'I'*400+b'\r\n')
            raw = b''.join(records).rstrip(b'\r\n')
            expected = self.split(raw, paired, 'python')
            self.assertEqual(self.split(raw, paired, 'native'), expected)
            self.assertEqual(self.split(raw, paired, 'native', replay=True), expected)

    def test_잘못된입력은같은유형으로거부(self):
        good=b'@r/1 HP:i:1\nACGT\n+\nIIII\n'
        cases=[(b'',False), (b'bad\nACGT\n+\nIIII\n',False),
               (b'@r\nACGT\n-\nIIII\n',False), (b'@r\nACGT\n+\nIII\n',False),
               (b'@r\nACGT\n+\n',False), (good,True),
               (good+good.replace(b'r/1',b'other/2'),True)]
        for raw,paired in cases:
            for mode in ['python','native']:
                with self.subTest(mode=mode,raw=raw),self.assertRaises(ValueError):
                    self.split(raw,paired,mode)

    def test_행상한_초과_거부(self):
        raw=b'@r\n'+b'A'*(64<<20)+b'\n+\nI\n'
        for mode in ['python','native']:
            with self.subTest(mode=mode),self.assertRaises(MemoryError):self.split(raw,False,mode)

    def test_낮은CPU_auto는Python_명시native는거부(self):
        with tempfile.TemporaryDirectory() as tmp:
            raw=b'@r\nA\n+\nI\n';paths=[Path(tmp)/'r.gz'];execution={}
            with patch.dict(os.environ,{'MFM_GPU_FASTQ_BACKEND':'auto'}):
                split_fastq(io.BytesIO(raw),paths,False,cpus=1,execution=execution)
            self.assertEqual(execution['backend'],'python')
            with patch.dict(os.environ,{'MFM_GPU_FASTQ_BACKEND':'native'}),self.assertRaises(ValueError):
                split_fastq(io.BytesIO(raw),paths,False,cpus=1)

    def test_auto는실측상불리한Long을Python으로유지(self):
        with tempfile.TemporaryDirectory() as tmp,patch.dict(os.environ,{'MFM_GPU_FASTQ_BACKEND':'auto'}):
            execution={}
            split_fastq(io.BytesIO(b'@r\nACGT\n+\nIIII\n'),[Path(tmp)/'r.gz'],False,cpus=68,execution=execution)
            self.assertEqual(execution['backend'],'python')

    def test_생산자예외는손실없이전파(self):
        class Broken:
            def read(self,size):raise OSError('입력 실패')
        with tempfile.TemporaryDirectory() as tmp,patch.dict(os.environ,{'MFM_GPU_FASTQ_BACKEND':'native'}):
            with self.assertRaisesRegex(OSError,'입력 실패'):
                split_fastq(Broken(),[Path(tmp)/'r.gz'],False,cpus=4)


if __name__=='__main__':unittest.main()
