from __future__ import annotations

import sys
import unittest
from unittest.mock import patch

from integrate_wgs_rna import cli


class Clair3QcExitTest(unittest.TestCase):
    """Stage 05가 QC JSON을 쓴 뒤에도 FAIL sample을 다음 stage에 넘기지 않는다."""

    @staticmethod
    def argv() -> list[str]:
        return [
            "mfm", "qc-longread-clair3", "--cohort", "T", "--output-root", "/tmp/out",
            "--platform", "hifi", "--model-path", "/tmp/model", "--reference-fasta", "/tmp/ref.fa",
            "--sample", "S1",
        ]

    def test_fail_or_missing_returns_nonzero(self) -> None:
        for result in ({"fail_count": 1, "missing_count": 0}, {"fail_count": 0, "missing_count": 1}):
            with patch.object(sys, "argv", self.argv()), patch.object(cli, "run_clair3_qc", return_value=result):
                with self.assertRaises(SystemExit) as raised:
                    cli.main()
                self.assertEqual(raised.exception.code, 2)

    def test_warn_only_remains_diagnostic(self) -> None:
        with patch.object(sys, "argv", self.argv()), patch.object(
            cli, "run_clair3_qc", return_value={"fail_count": 0, "missing_count": 0, "warn_count": 1}
        ):
            cli.main()


if __name__ == "__main__":
    unittest.main()
