"""동일 GRCh38 이름의 N/IUPAC 참조 혼용이 재발하지 않도록 검사한다."""
import hashlib
import gzip
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import pysam
from handoff import cli, bundled_assets, references, asset_restore


class HgsvcReferenceTests(unittest.TestCase):
    def test_project_callset_scope_and_missing_reference(self):
        refs={'autosomes':'ncbi.fa','hgsvc':'hgsvc.fa'}
        for project in ['HGSVC2','HGSVC3']:
            self.assertEqual(references.variant_reference({'project_id':project,'vcf_source':'public'},refs),'hgsvc.fa')
            with self.assertRaises(ValueError):
                references.variant_reference({'project_id':project,'vcf_source':'public'},{'autosomes':'ncbi.fa'})
        for project,source in [('HPRC_R2','public'),('1KGP_ILLUMINA','illumina_panel'),('HGSVC3','clair3')]:
            self.assertEqual(references.variant_reference({'project_id':project,'vcf_source':source},refs),'ncbi.fa')

    def test_bind_uses_hgsvc_for_sequence_and_vcf_but_retains_decode_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            config={'data_root':str(root),'state_root':str(root/'state'),'output_root':str(root/'out'),'run_id':'test'}
            row={'sample_id':'NA12878','platform':'hifi','project_id':'HGSVC2','vcf_source':'public',
                 'alignment_path':'wgs/HGSVC2/NA12878.bam','rna_path':'','data_group':'wgs','public_vcf_path':'public_vcf/HGSVC2/variants.vcf.gz'}
            refs={'ncbi':'decode.fa','autosomes':'ncbi.fa','hgsvc':'hgsvc.fa','broad':'broad.fa'}
            with patch.object(cli.planner,'read_tsv',return_value=[]),patch.object(cli.bundled_assets,'panel_directory',return_value=root/'panel'):
                path,bound=cli.bind_manifest(config,row,[],refs)
            self.assertTrue(path.is_file())
            self.assertEqual(bound['analysis_reference_fasta'],'hgsvc.fa')
            self.assertEqual(bound['vcf_reference_fasta'],'hgsvc.fa')
            self.assertEqual(bound['source_decode_reference_fasta'],'decode.fa')

    def test_bundled_hgsvc_is_required_only_for_hgsvc_and_validates_checksum(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);assets=root/'assets';config={'data_root':str(root/'data')}
            files=[]
            for name,content in [('GCA_000001405.15_GRCh38_no_alt_analysis_set.fna','>chr1\nC\n'),('hg38_autosomes.fa','>chr1\nC\n')]:
                p=assets/'references'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(content);pysam.faidx(str(p))
                for f in [p,Path(str(p)+'.fai')]:files.append({'path':str(f.relative_to(assets)),'size_bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
            manifest=assets/'manifest.json';manifest.write_text(json.dumps({'schema':'mfm_bundled_assets_v1','files':files}))
            self.assertIn('autosomes',bundled_assets.references(config,root))
            with self.assertRaises(ValueError):bundled_assets.references(config,root,needs_hgsvc=True)
            p=assets/'references/hgsvc_hg38_autosomes.fa';p.write_text('>chr1\nN\n');pysam.faidx(str(p))
            for f in [p,Path(str(p)+'.fai')]:files.append({'path':str(f.relative_to(assets)),'size_bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
            manifest.write_text(json.dumps({'schema':'mfm_bundled_assets_v1','files':files}))
            self.assertEqual(bundled_assets.references(config,root,needs_hgsvc=True)['hgsvc'],str(p))
            p.write_text('>chr1\nA\n')
            with self.assertRaises(ValueError):bundled_assets.references(config,root,needs_hgsvc=True)

    def test_strict_ref_check_passes_matching_n_and_rejects_wrong_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source.fa';wrong=root/'wrong.fa';source.write_text('>chr1\nACNTACGT\n');wrong.write_text('>chr1\nACGTACGT\n')
            for f in [source,wrong]:pysam.faidx(str(f))
            vcf=root/'callset.vcf';vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=8>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t3\t.\tN\tA\t.\tPASS\t.\n')
            for f,expected in [(source,0),(wrong,255)]:
                result=subprocess.run(['bcftools','norm','--check-ref','e','-f',str(f),'-Ou',str(vcf)],stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
                self.assertEqual(result.returncode,expected,result.stderr.decode())

    def test_unbundled_reference_uses_verified_source_and_keeps_generic_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);ref=root/'data/references';ref.mkdir(parents=True)
            ncbi=ref/'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
            ncbi.write_text(''.join('>chr'+str(i)+'\nACGT\n' for i in range(1,23)));pysam.faidx(str(ncbi))
            source=root/'source.fa.gz'
            with gzip.open(source,'wt') as handle:handle.write(''.join('>chr'+str(i)+'\nACNT\n' for i in range(1,23)))
            config={'data_root':str(root/'data'),'state_root':str(root/'state')}
            with patch.object(bundled_assets,'references',return_value=None),patch.object(references,'fetch',return_value=source) as fetch:
                refs=references.prepare(config,root,needs_hgsvc=True)
            args=fetch.call_args.args[0]
            self.assertEqual(args['checksum'],'6a5ceceed91ae5c14026454f09c4b6c4')
            self.assertEqual(args['size_bytes'],963047116)
            with pysam.FastaFile(refs['hgsvc']) as hgsvc,pysam.FastaFile(refs['autosomes']) as generic:
                self.assertEqual(hgsvc.fetch('chr1'),'ACNT');self.assertEqual(generic.fetch('chr1'),'ACGT')

    def test_lightweight_restore_preserves_n_and_reproduces_case_and_wrapping(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);expected=root/'expected.fa'
            expected.write_text(''.join(f'>chr{i}\n'+49*'A'+'N\n'+'NNNNN\n' for i in range(1,23)))
            pysam.faidx(str(expected));compressed=root/'source.gz'
            with gzip.open(compressed,'wt') as f:f.write(''.join(f'>chr{i}\n'+49*'a'+'nnnnnn\n' for i in range(1,23)))
            assets=root/'assets';records=[]
            for source,suffix in [(expected,''),(Path(str(expected)+'.fai'),'.fai')]:
                records.append({'path':asset_restore.HGSVC_AUTOSOMES+suffix,'size_bytes':source.stat().st_size,'sha256':hashlib.sha256(source.read_bytes()).hexdigest()})
            restorer=asset_restore.Restorer(root,assets,{'schema':'mfm_bundled_assets_v1','files':records},'samtools')
            with patch.object(asset_restore,'fetch',return_value=compressed) as fetch:
                result=restorer.restore(asset_restore.HGSVC_AUTOSOMES+'.fai')
            self.assertEqual(result.read_bytes(),Path(str(expected)+'.fai').read_bytes())
            self.assertEqual((assets/asset_restore.HGSVC_AUTOSOMES).read_bytes(),expected.read_bytes())
            self.assertEqual(fetch.call_args.args[0]['checksum'],'6a5ceceed91ae5c14026454f09c4b6c4')


if __name__=='__main__':unittest.main()
