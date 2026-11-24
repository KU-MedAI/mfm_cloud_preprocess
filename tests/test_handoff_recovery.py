"""실제 실패 유형인 curl56 재개, RNA 경로 격리, 참조 누락을 회귀 검증한다."""
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from handoff import download_retry, resume

class RecoveryTests(unittest.TestCase):
    def test_reset_reopens_range_from_updated_partial(self):
        with tempfile.TemporaryDirectory() as directory:
            partial = Path(directory) / 'input.part'
            offsets = []
            def run(command, **kwargs):
                offsets.append(partial.stat().st_size if partial.exists() else 0)
                with partial.open('ab') as handle:
                    handle.write(b'abc' if len(offsets) == 1 else b'def')
                return subprocess.CompletedProcess(command, 56 if len(offsets) == 1 else 0, '206')
            with patch.object(download_retry.subprocess, 'run', side_effect=run), patch.object(download_retry.time, 'sleep'):
                download_retry.transfer('ftp://example.test/input', partial,
                    lambda path: path.exists() and path.read_bytes() == b'abcdef')
            self.assertEqual(offsets, [0, 3])
            self.assertEqual(partial.read_bytes(), b'abcdef')

    def test_unknown_size_requires_successful_transfer(self):
        from handoff import downloads
        with tempfile.TemporaryDirectory() as directory:
            partial = Path(directory)/'input.bam.part'
            offsets = []
            def run(command, **kwargs):
                offsets.append(partial.stat().st_size if partial.exists() else 0)
                with partial.open('ab') as handle:handle.write(b'abc' if len(offsets)==1 else b'def')
                return subprocess.CompletedProcess(command,56 if len(offsets)==1 else 0,'206')
            resource = {'path':'input.bam','url':'ftp://example.test/input','size_bytes':None,'checksum':None}
            with patch.object(download_retry.subprocess,'run',side_effect=run),patch.object(download_retry.time,'sleep'):
                target = downloads.fetch(resource,directory,directory)
            self.assertEqual(offsets,[0,3])
            self.assertEqual(target.read_bytes(),b'abcdef')

    def test_http404_fails_without_retry(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(download_retry.subprocess, 'run', return_value=subprocess.CompletedProcess([], 22, '404')) as run:
                with self.assertRaises(subprocess.CalledProcessError):
                    download_retry.transfer('ftp://example.test/missing', Path(directory)/'part', lambda _:False)
                self.assertEqual(run.call_count, 1)

    def test_rna_replacement_preserves_url_and_original_row(self):
        with tempfile.TemporaryDirectory() as directory:
            row = {'sample_id':'HG00097', 'rna_path':'rna/old.bam'}
            config = {'data_root':directory, 'sample_rna_paths':{'HG00097':'recovery/HG00097.bam'}}
            original = [{'kind':'rna_alignment','path':'rna/old.bam','url':'https://example.test/original'},
                        {'kind':'rna_index','path':'rna/old.bam.bai','url':'https://example.test/original.bai'}]
            result = resume.sample_downloads(config,row,original)
            self.assertEqual(result[0]['path'], 'recovery/HG00097.bam')
            self.assertEqual(result[1]['path'], 'recovery/HG00097.bam.bai')
            self.assertEqual(result[0]['url'], original[0]['url'])
            self.assertEqual(original[0]['path'], 'rna/old.bam')
            config['sample_rna_paths']['HG00097'] = '../outside.bam'
            with self.assertRaises(ValueError):resume.rna_relative(config,row)

    def test_sample_resume_isolated_and_blacklist_required(self):
        with tempfile.TemporaryDirectory() as directory:
            config = {'pipeline_env':{'CLEAN_INTERMEDIATES_ON_SUCCESS':'0'},
                      'sample_pipeline_env':{'HG00126':{'WGS_START_AT':'11'}}}
            self.assertEqual(resume.sample_environment(config,'HG00126')['WGS_START_AT'], '11')
            self.assertNotIn('WGS_START_AT', resume.sample_environment(config,'HG00097'))
            rows = [{'sample_id':'HG00126'}]
            with self.assertRaises(ValueError):resume.require_blacklist(config,directory,rows)
            bed = Path(directory)/'references/hg38_autosomes_blacklist_union.bed'
            bed.parent.mkdir();bed.write_text('chr1\t0\t1\n')
            resume.require_blacklist(config,directory,rows)

if __name__ == '__main__':unittest.main()
