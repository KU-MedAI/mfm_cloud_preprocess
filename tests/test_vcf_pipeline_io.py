import json,os,signal,subprocess,sys,tempfile,unittest
from pathlib import Path
class VcfPipelineIoTests(unittest.TestCase):
 def run_child(self,pipeline):
  code="from integrate_wgs_rna.sample_vcf import run_pipeline;run_pipeline("+repr(pipeline)+",'IO_FIXTURE')"
  with subprocess.Popen([sys.executable,'-c',code],stdout=subprocess.PIPE,stderr=subprocess.PIPE,start_new_session=True) as child:
   try:out,err=child.communicate(timeout=5)
   except subprocess.TimeoutExpired:
    os.killpg(child.pid,signal.SIGKILL);child.communicate();raise
   return child.returncode,out,err
 def test_large_warning_does_not_block_output(self):
  with tempfile.TemporaryDirectory() as t:
   dest=Path(t)/'output';data='ACGT'*10000
   producer=[sys.executable,'-c',"import sys;sys.stderr.write('W'*200000);sys.stderr.flush();sys.stdout.write("+repr(data)+")"]
   consumer=[sys.executable,'-c',"import sys;open("+repr(str(dest))+",'wb').write(sys.stdin.buffer.read())"]
   code,out,err=self.run_child([producer,consumer]);self.assertEqual(code,0,err);self.assertEqual(dest.read_text(),data)
 def test_upstream_error_is_preserved_without_hanging(self):
  producer=[sys.executable,'-c',"import sys;sys.stderr.write('E'*200000+'ERROR_END');sys.exit(7)"]
  consumer=[sys.executable,'-c','import sys;sys.stdin.buffer.read()']
  code,out,err=self.run_child([producer,consumer]);self.assertNotEqual(code,0);self.assertIn(b'ERROR_END',err);self.assertIn(b'returncodes=[7]',err);self.assertGreater(len(err),200000)
 def test_middle_failure_and_downstream_failure_are_rejected(self):
  producer=[sys.executable,'-c',"import sys;sys.stdout.write('data')"]
  failure=[sys.executable,'-c',"import sys;sys.stderr.write('MIDDLE_END'+'E'*200000);sys.exit(9)"]
  consumer=[sys.executable,'-c','import sys;sys.stdin.buffer.read()']
  code,out,err=self.run_child([producer,failure,consumer]);self.assertNotEqual(code,0);self.assertIn(b'MIDDLE_END',err)
  code,out,err=self.run_child([producer,failure]);self.assertNotEqual(code,0);self.assertIn(b'MIDDLE_END',err)
 def test_spawn_failure_reaps_started_child(self):
  with tempfile.TemporaryDirectory() as t:
   code="""import subprocess,sys,time
from integrate_wgs_rna.sample_vcf import run_pipeline
started=[]
original=subprocess.Popen
def capture(*args,**kwargs):
 p=original(*args,**kwargs);started.append(p);return p
subprocess.Popen=capture
try:
 run_pipeline([[sys.executable,'-c','import time;time.sleep(60)'],['/missing/vcf-test-command']],'SPAWN_FIXTURE')
except FileNotFoundError:
 leaked=[p for p in started if p.poll() is None]
 for p in leaked:p.kill();p.wait()
 assert not leaked,'pipeline child leaked after spawn failure'
else:raise AssertionError('spawn failure was swallowed')
"""
   p=subprocess.run([sys.executable,'-c',code],capture_output=True,timeout=5);self.assertEqual(p.returncode,0,p.stderr)
if __name__=='__main__':unittest.main()
