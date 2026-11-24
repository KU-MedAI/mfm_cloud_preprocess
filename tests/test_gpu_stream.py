"""첫read·주석·paired순서·캐시·실패회수를 GPU 없이 전수 검증함."""
from pathlib import Path
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor
import gzip,json,tempfile,unittest
from integrate_wgs_rna.gpu_stream_inputs import StreamInputs
class StreamTests(unittest.TestCase):
 def setup_input(self,root,paired):
  source=root/'source.fq';parts=[[],[]] if paired else [[]]
  with source.open('wb') as f:
   for i in range(160):
    for mate in range(len(parts)):
     b=f'@r{i}/{mate+1} HP:i:1\tPS:i:17\tOD:i:1\n'.encode()+b'ACGT'*12500+b'\n+\n'+b'I'*50000+b'\n';parts[mate].append(b);f.write(b)
  script=root/'source.sh';script.write_text('set -e\ncat '+str(source)+'\n');return script,[b''.join(v) for v in parts]
 def test_첫record부터paired바이트와캐시일치함(self):
  for paired in [False,True]:
   with self.subTest(paired=paired),tempfile.TemporaryDirectory() as t:
    root=Path(t);src,expected=self.setup_input(root,paired);budget=SimpleNamespace(cpus=68,maximum_cpus=68,available_cpus=lambda:68,wait=lambda *a,**k:None,phase=lambda *a,**k:None)
    with StreamInputs(src,'key',root,paired,budget,3) as stream:
     self.assertFalse(stream.reused)
     with ThreadPoolExecutor(len(expected)) as pool:observed=list(pool.map(lambda p:p.read_bytes(),stream.paths))
     self.assertEqual(observed,expected);paths,count,reused=stream.finish();self.assertEqual(count,160*len(expected));self.assertFalse(reused)
     self.assertEqual([gzip.decompress(p.read_bytes()) for p in paths],expected)
    with StreamInputs(src,'key',root,paired,budget,3) as stream:self.assertTrue(stream.reused);self.assertEqual(stream.paths,paths)
 def test_source실패는완료cache를남기지않음(self):
  with tempfile.TemporaryDirectory() as t:
   root=Path(t);script=root/'bad.sh';script.write_text("printf '@bad\\nACGT\\n+\\nI\\n'\nexit 2\n");budget=SimpleNamespace(cpus=68,maximum_cpus=68,available_cpus=lambda:68,wait=lambda *a,**k:None,phase=lambda *a,**k:None)
   with self.assertRaises((ValueError,RuntimeError)):
    with StreamInputs(script,'bad',root,False,budget,3) as stream:
     stream.paths[0].read_bytes();stream.finish()
   self.assertFalse(list(root.glob('gpu_input_*/complete.json')))
 def test_소비자중단은producer를회수함(self):
  with tempfile.TemporaryDirectory() as t:
   root=Path(t);src,_=self.setup_input(root,False);budget=SimpleNamespace(cpus=68,maximum_cpus=68,available_cpus=lambda:68,wait=lambda *a,**k:None,phase=lambda *a,**k:None)
   with self.assertRaisesRegex(RuntimeError,'consumer abort'):
    with StreamInputs(src,'abort',root,False,budget,3) as stream:raise RuntimeError('consumer abort')
   self.assertFalse(stream.thread.is_alive());self.assertFalse(list(root.glob('gpu_input_*/complete.json')))
if __name__=='__main__':unittest.main()
