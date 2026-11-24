"""실제 FASTA/window 정책으로 서열·주석·mapping·재개·변조 검출을 검사함."""
import csv,json,os,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pysam,zarr
from integrate_wgs_rna.sequence_only import export_sequence,validate_contig,LAYOUT
from integrate_wgs_rna.window_index import WindowIndexConfig,run_window_index_generation
from integrate_wgs_rna.window_policy import validate_window_annotations


class SequenceOnlyTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name);self.root=self.base/'S1'
        self.sequence='NNN'+'ACGT'*25+'A'*45+'N'+'A'*54+'N'*100+'A'*12+'N'+'A'*12+'NNNN'
        for hap in ['HAP1','HAP2']:
            p=self.root/'individual_fasta'/hap/'chr1.fa';p.parent.mkdir(parents=True)
            p.write_text('>chr1\n'+self.sequence+'\n');pysam.faidx(str(p))
            n=len(self.sequence);p.with_suffix('.chain').write_text(f'chain 0 chr1 {n} + 0 {n} chr1 {n} + 0 {n} 1\n{n}\n')
        bed=self.base/'blacklist.bed';bed.write_text('chr1\t3\t103\n')
        run_window_index_generation(WindowIndexConfig(cohort='TEST',output_root=self.base,chroms=['chr1'],samples=['S1'],
             sequence_variant_policy='snv_indel',window_size=100,offset_unit=1,blacklist_bed=bed))
        self.directory=self.root/'compact'/LAYOUT/'chr1'

    def tearDown(self):self.temp.cleanup()

    def export(self):return export_sequence(self.root,contigs=['chr1'],jobs=1)

    def test_windowed_sequence_preserves_internal_n_blacklist_padding_and_chain(self):
        before=(self.root/'windows/windows.tsv').read_bytes()
        result=self.export();self.assertEqual(result['status'],'complete')
        self.assertEqual(before,(self.root/'windows/windows.tsv').read_bytes())
        self.assertFalse(any(self.directory.glob('readset*')));self.assertFalse((self.directory/'depth.zarr').exists())
        with (self.directory/'windows.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
        one=[r for r in rows if r['haplotype']=='HAP1']
        self.assertEqual([int(r['raw_start']) for r in one],[3,103,303])
        self.assertEqual([int(r['window_index']) for r in one],[0,1,2])
        self.assertEqual(one[0]['bl_frac'],'1.000000');self.assertEqual(one[1]['n_count'],'1');self.assertEqual(float(one[2]['n_fraction']),.04)
        g=zarr.open_group(str(self.directory/'hap_sequence.zarr'),mode='r');table=np.frombuffer(b'ACGTNP',dtype=np.uint8)
        for hap in ['HAP1','HAP2']:
            h=g[hap];offsets=h['offsets'][:]
            for i,row in enumerate(one):
                p=h['bases_packed'][offsets[i]:offsets[i+1]];codes=np.empty(p.size*2,dtype=np.uint8);codes[::2]=p>>4;codes[1::2]=p&15
                text=table[codes[:100]].tobytes().decode();start=int(row['raw_start']);actual=int(row['actual_length'])
                self.assertEqual(text,self.sequence[start:start+actual]+'P'*(100-actual))
        self.assertEqual(validate_contig(self.directory,source=True)['HAP1']['status'],'pass')
        self.assertEqual(validate_window_annotations(self.root)['status'],'pass')
        self.assertTrue((self.root/'metadata/reference_mapping/HAP1/chr1.chain.gz').is_file())
        trace=json.loads((self.root/'qc/stage12_trace/sequence_only/trace_manifest.json').read_text());self.assertEqual(trace['contigs'],{'chr1':'pass'})

    def test_resume_reuses_only_unchanged_input_and_rejects_replaced_fasta(self):
        self.export();chunk=self.directory/'hap_sequence.zarr/HAP1/bases_packed/0';stat=chunk.stat()
        self.assertTrue(self.export()['records']['chr1']['reused']);self.assertEqual(stat.st_mtime_ns,chunk.stat().st_mtime_ns)
        p=self.root/'individual_fasta/HAP1/chr1.fa';p.write_text(p.read_text().replace('ACGT','GCGT',1))
        with self.assertRaisesRegex(ValueError,'입력 지문'):self.export()

    def test_data_corruption_is_detected_before_resume_or_trace(self):
        self.export();a=zarr.open_array(str(self.directory/'hap_sequence.zarr/HAP1/bases_packed'),mode='a');a[0]=0x22
        with self.assertRaisesRegex(ValueError,'SHA|FASTA'):validate_contig(self.directory,source=True)
        with self.assertRaises(ValueError):self.export()

    def test_parallel_source_trace_does_not_publish_success_on_mismatch(self):
        from integrate_wgs_rna import sequence_only
        for hap in ('HAP1', 'HAP2'):
            fasta = self.root / 'individual_fasta' / hap / 'chr2.fa'
            fasta.write_text('>chr2\n' + self.sequence + '\n')
            pysam.faidx(str(fasta))
            chain = fasta.with_name('chr1.chain').read_text().replace('chr1', 'chr2')
            fasta.with_suffix('.chain').write_text(chain)
        contigs = ['chr1', 'chr2']
        run_window_index_generation(WindowIndexConfig(
            cohort='TEST', output_root=self.base, chroms=contigs, samples=['S1'],
            sequence_variant_policy='snv_indel', window_size=100, offset_unit=1,
            blacklist_bed=self.base / 'blacklist.bed',
        ))
        env = {'MFM_SAMPLE_CPU_BUDGET': '3', 'MFM_CPU_BUDGET': '3'}
        trace = self.root / 'qc/stage12_trace/sequence_only/trace_manifest.json'
        manifest = self.directory.parent / 'manifest.json'
        with patch.dict(os.environ, env):
            self.assertEqual(export_sequence(self.root, contigs=contigs)['status'], 'complete')
            previous_trace = trace.read_bytes()
            self.assertEqual(json.loads(previous_trace)['contigs'], {'chr1': 'pass', 'chr2': 'pass'})

            def alter_fasta_after_write(root):
                result = validate_window_annotations(root)
                fasta = root / 'individual_fasta/HAP1/chr2.fa'
                fasta.write_text(fasta.read_text().replace('ACGT', 'GCGT', 1))
                return result

            with patch.object(sequence_only, 'validate_window_annotations', alter_fasta_after_write):
                with self.assertRaisesRegex(ValueError, 'chr2 원본 FASTA 대조 실패'):
                    export_sequence(self.root, contigs=contigs)
            self.assertEqual(json.loads(manifest.read_text())['status'], 'running')
            self.assertEqual(trace.read_bytes(), previous_trace)

    def test_fake_read_or_depth_output_is_not_accepted(self):
        self.export();(self.directory/'depth.zarr').mkdir()
        with self.assertRaisesRegex(ValueError,'read/depth'):validate_contig(self.directory)

if __name__=='__main__':unittest.main()
