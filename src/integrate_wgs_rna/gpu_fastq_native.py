"""행별 Python 객체 생성을 피하는 bounded FASTQ 검증·분할 연결부."""
from threading import Thread
import subprocess
import tempfile

BATCH_BYTES = 4 << 20


def native_frames(tool, stream, paired, files, mirrors):
    errors = []
    with tempfile.TemporaryFile() as log:
        proc = subprocess.Popen([str(tool), str(int(paired))], stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=log)

        def feed():
            try:
                while True:
                    data = stream.read(BATCH_BYTES)
                    if not data:
                        break
                    proc.stdin.write(data)
            except BrokenPipeError:
                pass  # parser 오류는 아래 종료 코드/진단으로 보고한다.
            except BaseException as error:
                errors.append(error)
            finally:
                try:
                    proc.stdin.close()
                except BrokenPipeError:
                    pass

        thread = Thread(target=feed, daemon=True)
        thread.start()
        count = None
        try:
            while True:
                header = proc.stdout.read(5)
                if not header:
                    break
                if len(header) != 5 or count is not None:
                    raise RuntimeError('FASTQ 네이티브 프레임 손상')
                mate, size = header[0], int.from_bytes(header[1:], 'little')
                if not 0 < size <= BATCH_BYTES or (mate >= len(files) and not (mate == 2 and size == 8)):
                    raise RuntimeError('FASTQ 네이티브 프레임 범위 초과')
                data = proc.stdout.read(size)
                if len(data) != size:
                    raise RuntimeError('FASTQ 네이티브 프레임 잘림')
                if mate == 2:
                    count = int.from_bytes(data, 'little')
                else:
                    files[mate].write(data)
                    if mirrors:
                        mirrors[mate].write(data)
            code = proc.wait()
            thread.join(timeout=5)
            if errors:
                raise errors[0]
            log.seek(0)
            message = log.read(8192).decode(errors='replace').strip()
            if code == 75:
                raise MemoryError(message)
            if code == 65:
                raise ValueError(message)
            if code or count is None or thread.is_alive():
                raise RuntimeError('FASTQ 네이티브 검증 실패: ' + message)
            return count
        finally:
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    proc.kill(); proc.wait()
            proc.stdout.close()
            thread.join(timeout=1)
