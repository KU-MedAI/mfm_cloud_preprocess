from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from integrate_wgs_rna.full_preprocess import _replay_cli_commands


class ReplayCommandsTest(unittest.TestCase):
    def test_replay_command_uses_python_module_and_source_path(self) -> None:
        commands = _replay_cli_commands(
            ["/tmp/cli.py", "preprocess-full-shared-hap", "--sample-id", "A B"],
            python_executable="/env/jj_mfm/bin/python",
            module_root=Path("/project/src"),
        )
        self.assertEqual(commands[0], 'export PYTHONPATH=/project/src:"${PYTHONPATH:-}"')
        self.assertEqual(
            commands[1],
            "/env/jj_mfm/bin/python -m integrate_wgs_rna.cli preprocess-full-shared-hap --sample-id 'A B'",
        )


if __name__ == "__main__":
    unittest.main()
