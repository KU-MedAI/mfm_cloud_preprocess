from pathlib import Path
import importlib.util,sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from integrate_wgs_rna.model_input import dna_target_mask,masked_mean,torch_dna_cross_entropy
from integrate_wgs_rna.window_quality import depth_metrics

class ContractTest(unittest.TestCase):
    def test_n_padding_and_zero_loss(self):
        code=np.array([0,4,2,5,5],dtype=np.uint8)
        mask=dna_target_mask(code,3)
        self.assertEqual(mask.tolist(),[True,False,True,False,False])
        self.assertEqual(masked_mean([2,99,4,99,99],mask),3)
        self.assertEqual(masked_mean([1,2],[False,False]),0)
        self.assertEqual(depth_metrics([0,10,5,0,0],3),{'depth_sum':15,'depth_mean':5,'depth_zero_fraction':1/3,'depth_max':10})
        with self.assertRaises(ValueError):dna_target_mask(np.array([0,5,1]),3)

    @unittest.skipUnless(importlib.util.find_spec('torch'),'학습 환경의 optional torch adapter')
    def test_real_loss_gradient_ignores_n_and_padding(self):
        import torch
        logits=torch.zeros((1,5,4),requires_grad=True)
        codes=torch.tensor([[0,4,2,5,5]])
        loss=torch_dna_cross_entropy(logits,codes);loss.backward()
        self.assertTrue(torch.all(logits.grad[0,[1,3,4]]==0))
        self.assertTrue(torch.any(logits.grad[0,[0,2]]!=0))
        empty=torch.zeros((1,2,4),requires_grad=True)
        value=torch_dna_cross_entropy(empty,torch.tensor([[4,5]]));value.backward()
        self.assertEqual(value.item(),0);self.assertTrue(torch.all(empty.grad==0))


if __name__ == '__main__':unittest.main()
