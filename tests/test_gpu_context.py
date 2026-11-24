"""CuPy 초기화 장치 범위와 실패시 lease/mask 반환을 GPU 없이 검증함."""
from contextlib import ExitStack,nullcontext
from types import SimpleNamespace
from unittest.mock import patch
import os,unittest
from integrate_wgs_rna import gpu_arrays as m
class ContextTests(unittest.TestCase):
 def run_case(self,fail=False):
  closed=[];uuid='GPU-leased-2';lease=SimpleNamespace(device={'uuid':uuid,'pci':'0000:02:00.0'},reserved=4*(1<<30),close=lambda:closed.append(True))
  def devices():
   self.assertEqual(os.environ['CUDA_VISIBLE_DEVICES'],uuid)
   if fail:raise RuntimeError('CUDA init failed')
   return 1
  cp=SimpleNamespace(uint8='uint8',cuda=SimpleNamespace(runtime=SimpleNamespace(getDeviceCount=devices,deviceGetPCIBusId=lambda i:'0000:02:00.0'),Device=lambda i:nullcontext(),get_current_stream=lambda:SimpleNamespace(synchronize=lambda:None)),get_default_memory_pool=lambda:SimpleNamespace(set_limit=lambda **k:None),zeros=lambda *a,**k:SimpleNamespace(get=lambda:None))
  with ExitStack() as stack:
   stack.enter_context(patch.dict(os.environ,{'CUDA_VISIBLE_DEVICES':'0,1,2,3'}));stack.enter_context(patch.dict('sys.modules',{'cupy':cp}))
   for key,value in [('_lease',None),('_cp',None),('_owner',None),('_device',None),('_counts',{}),('_address_space',{})]:stack.enter_context(patch.object(m,key,value))
   stack.enter_context(patch('integrate_wgs_rna.gpu_pool.acquire',return_value=lease));stack.enter_context(patch.object(m.resource,'getrlimit',return_value=(-1,-1)))
   if fail:
    with self.assertRaisesRegex(RuntimeError,'CUDA init failed'):
     with m.session('zarr'):pass
    self.assertEqual(os.environ['CUDA_VISIBLE_DEVICES'],'0,1,2,3');self.assertTrue(closed)
   else:
    with m.session('zarr'):self.assertEqual(m.execution_metadata()['gpu_uuid'],uuid)
    m._close();self.assertTrue(closed)
 def test_예약GPU한장으로초기화함(self):self.run_case()
 def test_초기화실패는mask와예약을복원함(self):self.run_case(True)
if __name__=='__main__':unittest.main()
