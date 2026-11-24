"""실제 인덱스·STAR 버퍼·커널 상한 기반 RAM 예약을 검증한다."""
from pathlib import Path
import json
import os
import tempfile
import unittest
from unittest.mock import patch

from integrate_wgs_rna.rna_memory import command_memory_plan,estimate_star_memory,apply_star_memory
from integrate_wgs_rna.memory_policy import GIB,memory_limit,usable_memory


class RnaMemoryTests(unittest.TestCase):
    def test_실측인덱스와스레드버퍼를계산하고기존412GB보다작음(self):
        plan=estimate_star_memory('alignReads',68,index_bytes=29_627_155_261)
        self.assertEqual(plan['memory_bytes'],128*GIB)
        self.assertEqual(plan['components']['bam_sort_buffer'],29_627_155_261)
        self.assertGreaterEqual(plan['memory_bytes'],sum(plan['components'].values())*1.15)

    def test_커스텀큰sortRAM은예측에포함(self):
        plan=estimate_star_memory('alignReads',68,index_bytes=30_000_000_000,bam_sort_ram=400_000_000_000)
        self.assertGreater(plan['memory_bytes'],500_000_000_000)

    def test_index는생성상한과reference를반영(self):
        normal=estimate_star_memory('genomeGenerate',68,reference_bases=3_216_000_000)
        large=estimate_star_memory('genomeGenerate',68,reference_bases=9_000_000_000)
        self.assertEqual(normal['memory_bytes'],128*GIB)
        self.assertGreater(large['memory_bytes'],normal['memory_bytes'])

    def test_실제index파일크기를사용(self):
        with tempfile.TemporaryDirectory() as raw:
            directory=Path(raw)
            for name,size in [('Genome',100),('SA',800),('SAindex',50)]:
                with (directory/name).open('wb') as handle:handle.truncate(size)
            plan=command_memory_plan(['STAR','--runMode','alignReads','--genomeDir',raw,'--runThreadN','16'])
            self.assertEqual(plan['index_bytes'],950)
            self.assertEqual(plan['threads'],16)
            self.assertIsNone(command_memory_plan(['samtools','view','input.bam']))

    def test_동일leaseRAM만변경하고CPU와입력은보존(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();path=root/'leases/a.json'
            before={'branch':'RNA','stage':'16','sample':'S','cpus':32,'memory_bytes':384*GIB,'child_pid':123}
            path.write_text(json.dumps(before));plan=estimate_star_memory('alignReads',32,index_bytes=30_000_000_000)
            result=apply_star_memory(path,plan,running=True)
            self.assertEqual(result['memory_bytes'],128*GIB)
            self.assertEqual(result['cpus'],32)
            self.assertEqual(result['child_pid'],123)

    def test_실행중과소예약증가와전체예산초과는계산전거부(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();path=root/'leases/a.json'
            before={'branch':'RNA','cpus':32,'memory_bytes':64*GIB};path.write_text(json.dumps(before))
            plan=estimate_star_memory('alignReads',32,index_bytes=30_000_000_000)
            with self.assertRaisesRegex(ValueError,'실행 중'):apply_star_memory(path,plan,running=True)
            with patch('integrate_wgs_rna.rna_memory.usable_memory',return_value=64*GIB):
                with self.assertRaisesRegex(ValueError,'예산'):apply_star_memory(path,plan)
            self.assertEqual(json.loads(path.read_text()),before)

    def test_1점8TB기본상한과작은커널할당우선(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{},clear=True):
            root=Path(raw)
            with patch('integrate_wgs_rna.memory_policy.allocation_memory_limit',return_value=2_000_000_000_000):
                self.assertEqual(memory_limit(root),1_800_000_000_000)
                self.assertEqual(usable_memory(root),1_620_000_000_000)
            with patch('integrate_wgs_rna.memory_policy.allocation_memory_limit',return_value=1_717_986_918_400):
                self.assertEqual(memory_limit(root),1_717_986_918_400)
                self.assertEqual(usable_memory(root),1_546_188_226_560)
            with patch('integrate_wgs_rna.memory_policy.allocation_memory_limit',return_value=256*GIB):
                self.assertEqual(memory_limit(root),256*GIB)


if __name__=='__main__':unittest.main()
