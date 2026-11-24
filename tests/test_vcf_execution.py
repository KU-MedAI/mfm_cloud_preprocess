from pathlib import Path
import os,shutil,sys,tempfile,unittest
from unittest.mock import patch
import pysam
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from integrate_wgs_rna.sample_vcf import PerSampleVcfConfig,run_per_sample_vcf_extract


@unittest.skipUnless(shutil.which('bcftools'),'bcftools 필요')
class VcfExecutionTest(unittest.TestCase):
    def panels(self,root):
        root.mkdir()
        for chrom in ['chr1','chr2']:
            p=root/f'T.{chrom}.vcf'
            text='##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##contig=<ID=chr2,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n'
            for i,(ref,alt,gt) in enumerate([('A','C','0|1'),('A','C,G','0|2'),('AT','A','1|0'),('G','T','1|1'),('C','A','0|0'),('G','A','.|.')],1):
                text+=f'{chrom}\t{i}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\t{gt}\n'
            p.write_text(text);out=Path(str(p)+'.gz');pysam.tabix_compress(str(p),str(out),force=True);pysam.tabix_index(str(out),preset='vcf',force=True)

    def test_single_sequence_sample_releases_memory_only_after_worker_join(self):
        import json
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);self.panels(root/'panel');pool=root/'pool';(pool/'leases').mkdir(parents=True)
            lease=pool/'leases/stage.json'
            lease.write_text(json.dumps({'branch':'WGS','stage':'03','cpus':16,'memory_bytes':20<<30}))
            cfg=PerSampleVcfConfig(cohort='T',common_vcf_dir=root/'panel',output_dir=root/'out',
                 chroms=['chr1','chr2'],samples=['S1'],jobs=1,sequence_variant_policy='snv_indel')
            with patch.dict(os.environ,MFM_WGS_PAYLOAD='sequence_only',MFM_CPU_BUDGET='16',
                            MFM_CPU_POOL_LEASE_FILE=str(lease),MFM_RAM_LIMIT_BYTES=str(128<<30)):
                self.assertEqual(run_per_sample_vcf_extract(cfg)['ok_count'],1)
            row=json.loads(lease.read_text())
            self.assertEqual(row['cpus'],16)
            self.assertEqual(row['public_vcf_initial_memory_bytes'],20<<30)
            self.assertEqual(row['memory_bytes'],4<<30)
            self.assertEqual(row['public_vcf_memory_plan']['phase'],'sample_finalize')

    def test_new_mode_matches_legacy_and_cache_invalidates(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);self.panels(root/'panel')
            configs={mode:PerSampleVcfConfig(cohort='T',common_vcf_dir=root/'panel',output_dir=root/mode,
                      concat_dir=root/'concat',chroms=['chr1','chr2'],samples=['S1'],jobs=1,
                      sequence_variant_policy='snv_indel',extraction_mode=mode) for mode in ['chromosome','concat']}
            for mode in ['chromosome','concat']:run_per_sample_vcf_extract(configs[mode])
            for kind in ['seq','tag']:
                observed=[]
                for mode in ['chromosome','concat']:
                    with pysam.VariantFile(str(root/mode/f'S1.{kind}.vcf.gz')) as vcf:observed.append([str(x) for x in vcf])
                self.assertEqual(observed[0],observed[1])
            with patch('integrate_wgs_rna.sample_vcf.validate_tag_vcf',side_effect=AssertionError('캐시 미사용')):
                result=run_per_sample_vcf_extract(configs['chromosome'])
                self.assertEqual(result['skip_count'],1)
            self.assertFalse(list((root/'chromosome').glob('.S1.vcf.*')))
            # 여럿을 동시에 추출해도 같은 개인 변이 결과를 보존함.
            from dataclasses import replace
            multi=replace(configs['chromosome'],output_dir=root/'multi',samples=['S1','S2'],jobs=2)
            with patch.dict(os.environ,MFM_WGS_PAYLOAD='sequence_only'),patch('integrate_wgs_rna.public_vcf_parallel.reserve_stage_memory') as reserve:
                self.assertEqual(run_per_sample_vcf_extract(multi)['ok_count'],2)
                reserve.assert_not_called()
            with pysam.VariantFile(str(root/'multi/S2.seq.vcf.gz')) as second:
                self.assertGreater(len(list(second)),0)
            os.utime(root/'panel/T.chr1.vcf.gz',None)
            with self.assertRaisesRegex(RuntimeError,'입력/정책 변경'):
                run_per_sample_vcf_extract(configs['chromosome'])


if __name__ == '__main__':unittest.main()
