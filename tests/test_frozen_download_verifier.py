from __future__ import annotations

import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "verify_frozen", ROOT / "share_package/scripts/07_verify_frozen_downloads.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class FrozenDownloadVerifierTest(unittest.TestCase):
    def test_md5_and_multipart_etag(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "object.bin"
            path.write_bytes(b"abcdefghijkl")
            self.assertEqual(MODULE.md5_file(path), hashlib.md5(b"abcdefghijkl").hexdigest())
            parts = [hashlib.md5(chunk).digest() for chunk in (b"abcde", b"fghij", b"kl")]
            expected = f"{hashlib.md5(b''.join(parts)).hexdigest()}-3"
            self.assertEqual(MODULE.multipart_etag(path, part_size=5), expected)


if __name__ == "__main__":
    unittest.main()
