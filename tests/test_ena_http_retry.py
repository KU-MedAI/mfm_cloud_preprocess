"""ENA HTTPS403 복구가 범위·버전·부분 파일 보호를 유지하는지 검증함."""
from pathlib import Path
import io
import json
import sys
import tempfile
import unittest
from dataclasses import asdict
from unittest.mock import patch
from urllib.error import HTTPError
from urllib.request import Request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from handoff import http_download as http


ORIGINAL = 'https://ftp.sra.ebi.ac.uk/vol1/run/ERR1/read.bam'
MIRROR = ORIGINAL.replace('ftp.sra.ebi.ac.uk', 'ftp.ebi.ac.uk')


class Response(io.BytesIO):
    def __init__(self, data, content_range):
        super().__init__(data)
        self.status = 206
        self.headers = {'Content-Range': content_range, 'ETag': '"v1"',
                        'Last-Modified': 'Sat, 04 Mar 2023 14:00:45 GMT'}


class Tests(unittest.TestCase):
    def setUp(self):
        self.cache = patch.dict(http._EBI_HTTPS_MIRRORS, {}, clear=True)
        self.cache.start()
        self.addCleanup(self.cache.stop)

    def test_403시_HTTPS와_조건부헤더를유지하고_다음구간은바로대체호스트(self):
        request = Request(ORIGINAL, headers={'Range': 'bytes=4-7', 'If-Match': '"v1"',
                                           'Accept-Encoding': 'identity'})
        with patch.object(http, 'urlopen', side_effect=[HTTPError(ORIGINAL, 403, '', {}, None),
                                                       Response(b'efgh', 'bytes 4-7/8'),
                                                       Response(b'efgh', 'bytes 4-7/8')]) as opened:
            http._open(request, 12).close()
            http._open(request, 12).close()
        self.assertEqual([x.args[0].full_url for x in opened.call_args_list],
                         [ORIGINAL, MIRROR, MIRROR])
        self.assertEqual(opened.call_args_list[1].args[0].header_items(), request.header_items())
        self.assertEqual(opened.call_args_list[1].kwargs, {'timeout': 12})

    def test_다른HTTP오류와호스트와평문은자동전환하지않음(self):
        for url, code in [(ORIGINAL, 401), (ORIGINAL, 404), (ORIGINAL, 429), (ORIGINAL, 500),
                          (ORIGINAL.replace('https:', 'http:'), 403),
                          (ORIGINAL.replace('ftp.sra.ebi.ac.uk', 'other.example'), 403),
                          (ORIGINAL.replace('/vol1/', '/other/'), 403)]:
            with self.subTest(url=url, code=code):
                with patch.object(http, 'urlopen', side_effect=HTTPError(url, code, '', {}, None)) as opened:
                    with self.assertRaises(HTTPError):
                        http._open(Request(url), 12)
                    self.assertEqual(opened.call_count, 1)

    def test_대체실패는그대로보고하고_성공호스트로기억하지않음(self):
        with patch.object(http, 'urlopen', side_effect=[HTTPError(ORIGINAL, 403, '', {}, None),
                                                       HTTPError(MIRROR, 404, '', {}, None)]):
            with self.assertRaises(HTTPError):
                http._open(Request(ORIGINAL), 12)
        self.assertEqual(http._EBI_HTTPS_MIRRORS, {})

    def test_기존부분파일과원본URL장부를유지한채_이어받음(self):
        with tempfile.TemporaryDirectory() as directory:
            partial = Path(directory) / 'read.bam.part'
            partial.write_bytes(b'abcd')
            remote = http.RemoteObject(ORIGINAL, 8, '"v1"', 'Sat, 04 Mar 2023 14:00:45 GMT')
            metadata = Path(str(partial) + '.http.json')
            metadata.write_text(json.dumps(asdict(remote)))
            responses = [HTTPError(ORIGINAL, 403, '', {}, None),
                         Response(b'a', 'bytes 0-0/8'), Response(b'efgh', 'bytes 4-7/8')]
            with patch.object(http, 'urlopen', side_effect=responses), patch.dict('os.environ', {'MFM_HTTP_RANGE_WORKERS': '1'}):
                http.transfer(ORIGINAL, partial, lambda p: False, attempts=1, chunk_bytes=4)
            self.assertEqual(partial.read_bytes(), b'abcdefgh')
            self.assertEqual(json.loads(metadata.read_text()), asdict(remote))

    def test_대체주소에서도_버전변경을거부하고부분파일보존(self):
        with tempfile.TemporaryDirectory() as directory:
            partial = Path(directory) / 'read.bam.part'
            partial.write_bytes(b'abcd')
            remote = http.RemoteObject(ORIGINAL, 8, '"old"', None)
            metadata = Path(str(partial) + '.http.json')
            metadata.write_text(json.dumps(asdict(remote)))
            with patch.object(http, 'urlopen', side_effect=[HTTPError(ORIGINAL, 403, '', {}, None), Response(b'a', 'bytes 0-0/8')]):
                with self.assertRaises(http.RemoteChanged):
                    http.transfer(ORIGINAL, partial, lambda p: False, attempts=1)
            self.assertEqual(partial.read_bytes(), b'abcd')
            self.assertEqual(json.loads(metadata.read_text()), asdict(remote))


if __name__ == '__main__':
    unittest.main()
