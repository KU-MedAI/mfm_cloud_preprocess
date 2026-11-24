"""입력 경로와 별도로 출력 저장소를 검사하는 계약임."""
import tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from handoff.queue_resources import admit,output_space
class OutputSpaceTests(unittest.TestCase):
 def test_output_shortage_blocks_before_cpu_allocation(self):
  with tempfile.TemporaryDirectory() as t:
   c=dict(output_root=str(Path(t)/'future'),data_root='/tmp/input',min_output_free_bytes=100)
   with patch('handoff.queue_resources.os.statvfs',return_value=SimpleNamespace(f_bavail=80,f_frsize=1)) as stat:
    ok,reason=admit(c,t,[])
    self.assertFalse(ok);self.assertIn('최종 출력',reason);stat.assert_called_once_with(Path(t))
 def test_sufficient_output_space_passes(self):
  with patch('handoff.queue_resources.os.statvfs',return_value=SimpleNamespace(f_bavail=150,f_frsize=1)):
   self.assertTrue(output_space(dict(output_root='/tmp',min_output_free_bytes=100))[0])
 def test_stat_failure_blocks(self):
  with patch('handoff.queue_resources.os.statvfs',side_effect=OSError('unavailable')):
   self.assertFalse(output_space(dict(output_root='/tmp',min_output_free_bytes=100))[0])
 def test_unspecified_setting_preserves_existing_behavior(self):
  with patch('handoff.queue_resources.os.statvfs') as stat:
   self.assertTrue(output_space({})[0]);stat.assert_not_called()
if __name__=='__main__':unittest.main()
