"""실제 bcftools의 public VCF 직렬/병렬 레코드·GT·PS·실패 계약을 검사한다."""
from pathlib import Path
from dataclasses import replace
import contextlib
import io
import json
import os
import shutil
import tempfile
import unittest
from unittest.mock import patch
import pysam
from integrate_wgs_rna.public_vcf import PublicVcfConfig, prepare_public_vcf
from integrate_wgs_rna.public_vcf_parallel import resource_plan, reserve_stage_memory


class PublicVcfParallelTests(unittest.TestCase):
    def fixture(self,raw):
        root=Path(raw);ref=root/'ref.fa';ref.write_text('>chr1\n'+'A'*1000+'\n>chr2\n'+'A'*1000+'\n');pysam.faidx(str(ref))
        header='##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n##FORMAT=<ID=PS,Number=1,Type=Integer,Description="PS">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSOURCE\n'
        sources=[]
        for index in range(2):
            p=root/f'source{index}.vcf'
            rows=[]
            for chrom in ['chr1','chr2']:
                rows += [f'{chrom}\t10\tr{index}\tA\tC,G\t30\tPASS\t.\tGT:PS\t1|2:{100+index}\n',
                         f'{chrom}\t20\ti{index}\tAA\tA\t40\tPASS\t.\tGT:PS\t0|1:{200+index}\n',
                         f'{chrom}\t30\th{index}\tA\tT\t20\tPASS\t.\tGT:PS\t1|1:.\n']
            p.write_text(header+''.join(rows));gz=Path(str(p)+'.gz');pysam.tabix_compress(str(p),str(gz));pysam.tabix_index(str(gz),preset='vcf');sources.append(gz)
        config=PublicVcfConfig(sample='ALIAS',source_sample='SOURCE',source_vcfs=sources,reference_fasta=ref,
            output_dir=root/'serial',callset='fixture',phase_source='provided',chroms=['chr1','chr2'],bcftools_bin=shutil.which('bcftools'))
        return config

    def records(self,path):
        with pysam.VariantFile(str(path)) as vcf:return [str(record) for record in vcf]

    def test_직렬병렬의정렬_GT_PS_중복_INDEL_alias통계동일(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{'MFM_CPU_BUDGET':'8','MFM_CPU_POOL_GRANTED':'8'}):
            config=self.fixture(raw)
            with patch.dict(os.environ,{'MFM_PUBLIC_VCF_JOBS':'1'}):first=prepare_public_vcf(config)
            parallel=replace(config,output_dir=Path(raw)/'parallel')
            stream=io.StringIO()
            with contextlib.redirect_stdout(stream),patch.dict(os.environ,{'MFM_PUBLIC_VCF_JOBS':'auto'}):second=prepare_public_vcf(parallel)
            for name in ['ALIAS.seq.vcf.gz','ALIAS.tag.vcf.gz']:
                self.assertEqual(self.records(config.output_dir/name),self.records(parallel.output_dir/name))
            for key in ['seq_stats','tag_stats']:
                if key in first:self.assertEqual(first[key],second[key])
            completed=[json.loads(line.split('] ',1)[1]) for line in stream.getvalue().splitlines() if line.startswith('[PUBLIC-VCF-PARALLEL]')][-1]
            self.assertEqual(completed['completed_chromosomes'],2)
            self.assertEqual(completed['peak_active_chromosomes'],2)
            self.assertFalse(list(parallel.output_dir.glob('.public_vcf_parallel_*')))

    def test_REF오류는실패하고부분결과와임시샤드를남기지않음(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{'MFM_CPU_BUDGET':'8','MFM_CPU_POOL_GRANTED':'8'}):
            config=self.fixture(raw);config.reference_fasta.write_text('>chr1\n'+'C'*1000+'\n>chr2\n'+'C'*1000+'\n');pysam.faidx(str(config.reference_fasta))
            with self.assertRaises(Exception):prepare_public_vcf(config)
            self.assertFalse(list(config.output_dir.glob('*.vcf.gz')))
            self.assertFalse(list(config.output_dir.glob('.public_vcf_parallel_*')))

    def test_CPU와RAM에서동시염색체수를정하고명시직렬존중(self):
        with patch.dict(os.environ,{'MFM_CPU_BUDGET':'68','MFM_CPU_POOL_GRANTED':'8','MFM_PUBLIC_VCF_JOBS':'auto'}):
            plan=resource_plan(22);self.assertLessEqual(plan['workers']*2+1,8)
            self.assertLessEqual(plan['estimated_peak_bytes'],plan['ram_budget_bytes'])
            with patch.dict(os.environ,{'MFM_PUBLIC_VCF_JOBS':'1'}):self.assertEqual(resource_plan(22)['workers'],1)

    def test_단계피크로감축해도CPU보존하며예산초과증액거부(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();path=root/'leases/self.json'
            row={'branch':'WGS','stage':'03','cpus':27,'memory_bytes':128*1024**3}
            path.write_text(json.dumps(row))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path),'MFM_CPU_BUDGET':'68','MFM_CPU_POOL_GRANTED':'27'}):
                plan=resource_plan(22);reserve_stage_memory(plan)
                self.assertEqual(resource_plan(22)['workers'],plan['workers'])
            actual=json.loads(path.read_text());self.assertEqual(actual['cpus'],27)
            self.assertEqual(actual['public_vcf_initial_memory_bytes'],128*1024**3)
            self.assertEqual(actual['memory_bytes'],plan['estimated_peak_bytes'])
            with self.assertRaises(ValueError):reserve_stage_memory(dict(plan,estimated_peak_bytes=129*1024**3),path)


    def test_태그생성전_CPU_RAM반환과캐시재사용(self):
        from integrate_wgs_rna.public_vcf import build_tag_vcf
        from integrate_wgs_rna.public_vcf_parallel import entry_memory
        with tempfile.TemporaryDirectory() as raw:
            config=self.fixture(raw);root=Path(raw)/'pool';(root/'leases').mkdir(parents=True)
            lease=root/'leases/self.json'
            lease.write_text(json.dumps({'branch':'WGS','stage':'03','cpus':8,'memory_bytes':entry_memory(8)}))
            observed=[]
            def tag(*args):
                observed.append(json.loads(lease.read_text()))
                return build_tag_vcf(*args)
            env={'MFM_CPU_BUDGET':'8','MFM_CPU_POOL_GRANTED':'8','MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_PUBLIC_VCF_JOBS':'auto'}
            with patch.dict(os.environ,env),patch('integrate_wgs_rna.public_vcf.build_tag_vcf',side_effect=tag):
                first=prepare_public_vcf(config)
                cached=prepare_public_vcf(config)
            self.assertEqual(len(observed),1)
            self.assertEqual(observed[0]['cpus'],2)
            self.assertEqual(observed[0]['memory_bytes'],4*1024**3)
            self.assertEqual(observed[0]['public_vcf_memory_plan']['phase'],'tag_and_validation')
            self.assertEqual(cached['action'],'skip')
            self.assertEqual(first['seq_stats'],cached['seq_stats'])
            self.assertEqual(first['tag_stats'],cached['tag_stats'])

    def test_검증예약누락은조용히진행하지않음(self):
        from integrate_wgs_rna.public_vcf_parallel import reserve_validation_memory
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(Path(raw)/'missing.json')}):
            with self.assertRaises(FileNotFoundError):reserve_validation_memory()



if __name__=='__main__':unittest.main()
