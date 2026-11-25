"""공통36채널·HP 분리·BAF 분모·native/Python/HTSlib 경계를 검증함."""
from pathlib import Path
import hashlib
import os
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT/'src'), str(ROOT/'vendor/MFM-DNARNAFM')]
from integrate_wgs_rna.pileup_features import compute as legacy
from integrate_wgs_rna.pileup_features.export import filter_config
from integrate_wgs_rna.pileup_features.complete_compute import ReferenceCache, calculate_counts, projected
from integrate_wgs_rna.pileup_features.complete_native import CompleteNative
from integrate_wgs_rna.pileup_features.projection import Projection
from integrate_wgs_rna.pileup_features.schema import CHANNELS, LEGACY_CHANNELS
from integrate_wgs_rna.pileup_features.variants import VariantLookup

CI = {name: i for i, name in enumerate(CHANNELS)}


def write_bam(path, reference, records):
    header = {'HD': {'VN':'1.6', 'SO':'coordinate'},
              'SQ':[{'SN':'chr1', 'LN':len(reference), 'M5':hashlib.md5(reference.encode()).hexdigest()}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for i, record in enumerate(sorted(records, key=lambda r:r.get('start',0))):
            read = pysam.AlignedSegment()
            read.query_name = 'r'+str(i)
            read.query_sequence = record['seq']
            read.reference_id = 0
            read.reference_start = record.get('start',0)
            read.flag = record.get('flag',0)
            read.mapping_quality = record.get('mapq',60)
            read.cigartuples = record.get('cigar',[(0,len(record['seq']))])
            if record.get('bq',30) is not None:
                read.query_qualities = [record.get('bq',30)]*len(record['seq'])
            if 'hp' in record:
                read.set_tag('HP',record['hp'])
            out.write(read)
    pysam.index(str(path))


class CompletePileupTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.sequence = 'ACGT'*16
        self.ref = self.root/'ref.fa'
        self.ref.write_text('>chr1\n'+self.sequence+'\n')
        pysam.faidx(str(self.ref))
        self.vcf = self.make_vcf(['chr1\t1\t.\tA\tG\t60\tPASS\t.\tGT\t0|1'])

    def tearDown(self):
        self.temp.cleanup()

    def make_vcf(self, records):
        path = self.root/'variants.vcf'
        path.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=64>\n'
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
                        +'\n'.join(records)+'\n')
        gz = Path(str(path)+'.gz')
        pysam.tabix_compress(str(path),str(gz),force=True)
        pysam.tabix_index(str(gz),preset='vcf',force=True)
        return gz

    def values(self, records, *, read_type='short_read', oracle=False, block=17):
        path = self.root/'reads.bam'
        write_bam(path,self.sequence,records)
        with pysam.AlignmentFile(str(path),'rb') as bam:
            cache = ReferenceCache(bam,'chr1',filter_config(read_type,30),self.ref,self.vcf,'S',oracle=oracle,block=block)
            try:
                return cache.gather(np.arange(len(self.sequence)))
            finally:
                cache.close()

    def test_schema_places_reverse_after_forward_and_validates_both_orders(self):
        from integrate_wgs_rna.pileup_features.schema import FEATURE_SET, PREVIOUS_FEATURE_SET, PREVIOUS_CHANNELS, resolve_channels
        expected = ('log1p_depth','mean_mapq','mean_baseq','strand_forward_ratio','strand_reverse_ratio',
                    'read_A_ratio','read_C_ratio','read_G_ratio','read_T_ratio','read_N_ratio',
                    'insertion_ratio','deletion_after_event_ratio','deletion_after_max_len',
                    'hp1_log1p_depth','hp1_read_A_ratio','hp1_read_C_ratio','hp1_read_G_ratio','hp1_read_T_ratio',
                    'hp1_read_N_ratio','hp1_insertion_ratio','hp1_deletion_after_event_ratio','hp1_deletion_after_max_len',
                    'hp2_log1p_depth','hp2_read_A_ratio','hp2_read_C_ratio','hp2_read_G_ratio','hp2_read_T_ratio',
                    'hp2_read_N_ratio','hp2_insertion_ratio','hp2_deletion_after_event_ratio','hp2_deletion_after_max_len',
                    'vcf_ref_allele_fraction','baf','minor_allele_fraction',
                    'read_supports_variant_ratio','base_diff_from_ref_ratio')
        self.assertEqual(CHANNELS,expected)
        self.assertEqual(len(CHANNELS),36)
        self.assertEqual(len(set(CHANNELS)),36)
        self.assertEqual(tuple(name for name in CHANNELS if name in LEGACY_CHANNELS),LEGACY_CHANNELS)
        self.assertEqual(PREVIOUS_CHANNELS[:12],LEGACY_CHANNELS)
        self.assertEqual(PREVIOUS_CHANNELS[35],'strand_reverse_ratio')
        for channels, feature_set in [(CHANNELS,FEATURE_SET),(PREVIOUS_CHANNELS,PREVIOUS_FEATURE_SET)]:
            self.assertEqual(resolve_channels(list(channels),feature_set),channels)
        for channels, feature_set in [(CHANNELS,PREVIOUS_FEATURE_SET),(PREVIOUS_CHANNELS,FEATURE_SET),(CHANNELS,'unknown')]:
            with self.assertRaisesRegex(ValueError,'채널 순서'):resolve_channels(list(channels),feature_set)
        self.assertEqual(CHANNELS.count('strand_forward_ratio'),1)
        self.assertIn('strand_reverse_ratio',CHANNELS)
        self.assertNotIn('hp1_mean_baseq',CHANNELS)
        self.assertNotIn('hp1_baf',CHANNELS)
        self.assertNotIn('hp1_strand_forward_ratio',CHANNELS)

    def test_new_execution_defaults_to_pileup_and_explicit_readset_is_preserved(self):
        from handoff.sequence_mode import output_mode,runtime_environment
        from integrate_wgs_rna.pileup_features.schema import LAYOUT
        for read_type in ('short_read','long_read'):
            row=dict(sample_id='S',read_type=read_type,rna_path='')
            self.assertEqual(output_mode({},row),'pileup_feature')
            self.assertEqual(runtime_environment({},row)['COMPACT_LAYOUT_NAME'],LAYOUT)
            config=dict(pipeline_env={'MFM_SHORT_READ_OUTPUT_MODE':'full','MFM_LONG_READ_OUTPUT_MODE':'full'})
            self.assertEqual(output_mode(config,row),'full')

    def test_hp_evidence_and_baf_are_different_denominators(self):
        records = [dict(seq='ACGT',hp=1)]*6 + [dict(seq='GCGT',hp=2,flag=16)]*4
        records += [dict(seq='NCGT',bq=None)]*2 + [dict(seq='TCGT',flag=16)]
        values = self.values(records)
        f = values['features'][0]
        self.assertEqual(int(values['depth'][0]),13)
        self.assertEqual(int(values['hp1_depth'][0]),6)
        self.assertEqual(int(values['hp2_depth'][0]),4)
        for key, expected in [('read_A_ratio',6/13),('read_G_ratio',4/13),('read_N_ratio',2/13),
                              ('hp1_read_A_ratio',1),('hp2_read_G_ratio',1),
                              ('vcf_ref_allele_fraction',.6),('baf',.4),('minor_allele_fraction',.4),
                              ('read_supports_variant_ratio',4/13),('base_diff_from_ref_ratio',7/13),
                              ('strand_forward_ratio',8/13),('strand_reverse_ratio',5/13)]:
            self.assertAlmostEqual(float(f[CI[key]]),expected,places=6,msg=key)
        self.assertEqual(int(values['baf_valid'][0]),1)
        self.assertEqual(int(values['baf_heterozygous'][0]),1)
        self.assertEqual(int(values['baseq_depth'][0]),11)
        self.assertTrue(np.all(values['features'][1:,CI['baf']]==0))

    def test_native_python_and_independent_htslib_match_complex_cigars(self):
        records = [dict(seq='ACTGT',cigar=[(0,2),(1,1),(0,2)],hp=1),
                   dict(seq='ACGT',cigar=[(0,2),(2,2),(0,2)],hp=2,flag=16,bq=None),
                   dict(seq='TACGT',cigar=[(4,1),(0,2),(6,1),(0,2)],hp=0),
                   dict(seq='ACGT',cigar=[(0,2),(3,2),(0,2)],hp='1')]
        records += [dict(seq='ACGT',flag=flag,hp=1,mapq=mq) for flag,mq in
                    [(0,29),(256,60),(512,60),(1024,60),(2048,60)]]
        path = self.root/'complex.bam';write_bam(path,self.sequence,records)
        for read_type in ('short_read','long_read'):
            with pysam.AlignmentFile(str(path),'rb') as bam:
                config = filter_config(read_type,30)
                expected = calculate_counts(bam,'chr1',0,64,config)
                for oracle in (False,True):
                    native = CompleteNative(bam,'chr1',config,oracle=oracle)
                    try:
                        actual = native.calculate(0,64)
                        for a,b in zip(actual,expected):np.testing.assert_array_equal(a,b)
                        left,right = native.calculate(0,3),native.calculate(3,64)
                        for i in (0,1):np.testing.assert_array_equal(actual[i],np.concatenate([left[i],right[i]],axis=2))
                    finally:native.close()

    def test_first_twelve_match_legacy_native_and_empty_baf_has_mask(self):
        records = [dict(seq='ACGT',hp=1),dict(seq='AC',cigar=[(0,2),(2,2)],hp=2),
                   dict(seq='GCGT',bq=None),dict(seq='ACGT',flag=2048,hp=1)]
        values = self.values(records,read_type='long_read')
        path = self.root/'reads.bam'
        with pysam.AlignmentFile(str(path),'rb') as bam:
            old = legacy.ReferenceCache(bam,'chr1',filter_config('long_read',30))
            try:expected = old.gather(np.arange(64))
            finally:old.close()
        np.testing.assert_array_equal(values['features'][:,[CI[name] for name in LEGACY_CHANNELS]],expected['features'])
        for key in ('depth','baseq_depth','supplementary_depth','max_deletion_length'):
            np.testing.assert_array_equal(values[key],expected[key])
        self.assertEqual(int(values['baf_valid'][20]),0)

    def test_variant_ambiguity_missing_gt_and_reference_mismatch(self):
        self.vcf = self.make_vcf(['chr1\t1\t.\tA\tG\t60\tPASS\t.\tGT\t0|1',
                                 'chr1\t1\t.\tA\tT\t60\tPASS\t.\tGT\t0|1',
                                 'chr1\t2\t.\tC\tT,G\t60\tPASS\t.\tGT\t1|2',
                                 'chr1\t3\t.\tG\tA\t60\tPASS\t.\tGT\t.|.'])
        values = self.values([dict(seq='ACGT',hp=1)])
        self.assertFalse(values['baf_site'].any())
        with self.assertRaisesRegex(ValueError,'sample'):
            VariantLookup(self.vcf,'WRONG')
        self.vcf = self.make_vcf(['chr1\t1\t.\tC\tG\t60\tPASS\t.\tGT\t0|1'])
        with self.assertRaisesRegex(ValueError,'REF/reference'):
            self.values([dict(seq='ACGT')])

    def test_insertion_tokens_keep_event_support_but_no_fake_baf_or_base_depth(self):
        path = self.root/'insert.bam'
        write_bam(path,self.sequence,[dict(seq=self.sequence[:4]+'TT'+self.sequence[4:],cigar=[(0,4),(1,2),(0,60)],hp=1)])
        chain = self.root/'hap.chain'
        chain.write_text('chain 1 chr1 64 + 0 64 chr1 66 + 0 66 1\n4 0 2\n60\n')
        projection = Projection(chain,'chr1',64,66)
        with pysam.AlignmentFile(str(path),'rb') as bam:
            cache = ReferenceCache(bam,'chr1',filter_config('short_read',30),self.ref,self.vcf,'S',block=8)
            try:values=projected(cache,projection,0,66,70)
            finally:cache.close()
        np.testing.assert_array_equal(values['mapping_kind'][4:6],[2,2])
        self.assertTrue(np.all(values['features'][4:6,CI['insertion_ratio']]==1))
        self.assertTrue(np.all(values['features'][4:6,CI['hp1_insertion_ratio']]==1))
        for name in ('depth','hp1_depth','hp2_depth','baf_valid','baf_site'):
            self.assertFalse(values[name][4:6].any(),name)
        self.assertTrue(np.all(values['reference_base_code'][4:6]==255))
        self.assertFalse(values['features'][66:].any())


class CompletePackageTest(unittest.TestCase):
    def setUp(self):
        from test_pileup_feature import PileupPackageTest
        from integrate_wgs_rna.pileup_features.schema import LAYOUT
        PileupPackageTest.setUp(self)
        self.directory=self.root/'compact'/LAYOUT/'chr1'
        path=self.base/'S1.seq.vcf'
        path.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length='+str(len(self.sequence))+'>\n'
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
                        'chr1\t3\t.\tA\tG\t60\tPASS\t.\tGT\t0|0\n')
        self.vcf=Path(str(path)+'.gz')
        pysam.tabix_compress(str(path),str(self.vcf),force=True)
        pysam.tabix_index(str(self.vcf),preset='vcf',force=True)

    def tearDown(self):self.temp.cleanup()

    def run_export(self):
        from integrate_wgs_rna.pileup_features.export import export
        with patch.dict(os.environ,{'MFM_SAMPLE_CPU_BUDGET':'1','MFM_RNA_OUTPUT_MODE':'coverage_only'}):
            return export(self.root,self.bam,self.ref,self.ref,contigs=['chr1'],read_type='short_read',
                          variants=self.vcf,variant_sample='S1',jobs=1)

    def test_export_reader_resume_and_cleanup_preserve_36_channel_contract(self):
        import json,zarr
        from integrate_wgs_rna.pileup_features.export import validate_contig
        from integrate_wgs_rna.pileup_features.schema import FEATURE_SET
        from integrate_wgs_rna.window_reader import WindowReader
        from handoff.cleanup import cleanup
        before=(self.root/'windows/windows.tsv').read_bytes()
        report=self.run_export()
        self.assertEqual(report['feature_set'],FEATURE_SET)
        self.assertEqual(before,(self.root/'windows/windows.tsv').read_bytes())
        group=zarr.open_group(str(self.directory/'pileup.zarr'),mode='r')
        self.assertEqual(group['HAP1/features'].shape[1],36)
        self.assertEqual(group['HAP1/features'].order,'F')
        self.assertEqual(group['HAP1/features'].chunks,(16,36))
        self.assertEqual(group['HAP1/depth'].chunks,(16,))
        self.assertFalse(group['HAP1/baf_heterozygous'][:].any())
        self.assertEqual(group['HAP1/baf_heterozygous'].nchunks_initialized,0)
        self.assertTrue(group['HAP1/hp1_depth'][:].any())
        self.assertTrue(group['HAP1/hp2_depth'][:].any())
        path=self.directory/'pileup.zarr/HAP1/features/0.0';mtime=path.stat().st_mtime_ns
        self.run_export();self.assertEqual(mtime,path.stat().st_mtime_ns)
        from test_pileup_feature import PileupPackageTest
        PileupPackageTest._finalize(self)
        report=json.loads((self.root/'qc/final_package_validation.json').read_text())
        self.assertEqual(report['pileup_feature_set'],FEATURE_SET)
        reader=WindowReader(self.root)
        try:
            value=reader.load('chr1','HAP2',0)
            self.assertEqual(value['pileup_features'].shape,(16,36))
            self.assertTrue(value['pileup_features'].flags.c_contiguous)
            self.assertEqual(tuple(value['pileup_channels']),CHANNELS)
            self.assertEqual(value['pileup_channels'][4],'strand_reverse_ratio')
        finally:reader.close()
        cleanup(self.root,'S1','wgs',output_root=self.base)
        self.assertFalse((self.root/'individual_fasta').exists())
        self.assertEqual(validate_contig(self.directory)['HAP1']['status'],'pass')
        reader=WindowReader(self.root)
        try:self.assertEqual(reader.load('chr1','HAP1',0)['pileup_features'].shape[1],36)
        finally:reader.close()

    def test_previous_36_order_reads_and_validates_without_relabel_or_overwrite(self):
        import json,time,zarr
        from integrate_wgs_rna.pileup_features.schema import PREVIOUS_CHANNELS, PREVIOUS_FEATURE_SET, POLICY
        from integrate_wgs_rna.pileup_features.export import validate_contig
        from integrate_wgs_rna.pipeline.checkpoints import digest
        from integrate_wgs_rna.window_reader import WindowReader
        from handoff.sequence_mode import payload_matches
        self.run_export()
        group=zarr.open_group(str(self.directory/'pileup.zarr'),mode='a')
        order=[CI[name] for name in PREVIOUS_CHANNELS]
        original=group['HAP2/features'][:16].copy()
        # 이전 발행본의 열·정책·지문을 가진 작은 패키지를 재현함.
        path=self.directory/'contig_metadata.json';metadata=json.loads(path.read_text())
        policy=dict(POLICY,feature_set=PREVIOUS_FEATURE_SET)
        provenance=metadata['pileup']['provenance']
        provenance.update(channels=list(PREVIOUS_CHANNELS),producer=PREVIOUS_FEATURE_SET,policy=policy)
        fingerprint=digest(provenance)
        metadata['pileup']['fingerprint_sha256']=fingerprint
        for hap in ('HAP1','HAP2'):
            # 이전 C-order·작은 chunk 저장도 그대로 읽고 source 검증해야 함.
            previous=group[hap]['features']
            values=previous[:][:,order];codec=previous.compressor
            del group[hap]['features']
            array=group[hap].create_dataset('features',data=values,chunks=(8,36),order='C',compressor=codec)
            metadata['pileup']['content_sha256'][hap]['features']=hashlib.sha256(array[:].tobytes()).hexdigest()
        path.write_text(json.dumps(metadata))
        group.attrs.update(channels=list(PREVIOUS_CHANNELS),feature_set=PREVIOUS_FEATURE_SET,policy=policy,fingerprint_sha256=fingerprint)
        manifest=self.directory.parent/'manifest.json'
        for path in (manifest,self.root/'metadata/preprocess.json'):
            data=json.loads(path.read_text());data.update(channels=list(PREVIOUS_CHANNELS),feature_set=PREVIOUS_FEATURE_SET,feature_policy=policy)
            path.write_text(json.dumps(data))
        path=self.root/'metadata/sample.json';data=json.loads(path.read_text())
        data.update(pileup_feature_set=PREVIOUS_FEATURE_SET,pileup_feature_channels=list(PREVIOUS_CHANNELS));path.write_text(json.dumps(data))
        metrics={}
        validate_contig(self.directory,source=True,metrics=metrics)
        path=self.root/'qc/stage12_trace/pileup_feature/trace_manifest.json'
        trace=json.loads(path.read_text());trace.update(source_validation={'chr1':metrics},completed_at=time.time())
        path.write_text(json.dumps(trace))
        from test_pileup_feature import PileupPackageTest
        PileupPackageTest._finalize(self)
        report=json.loads((self.root/'qc/final_package_validation.json').read_text())
        self.assertEqual(report['pileup_feature_set'],PREVIOUS_FEATURE_SET)
        reader=WindowReader(self.root)
        try:
            value=reader.load('chr1','HAP2',0)
            self.assertEqual(tuple(value['pileup_channels']),PREVIOUS_CHANNELS)
            np.testing.assert_array_equal(value['pileup_features'],original[:,order])
        finally:reader.close()
        row=dict(sample_id='S1',read_type='short_read',rna_path='')
        self.assertFalse(payload_matches({},row,self.root))
        before=manifest.read_bytes()
        with self.assertRaisesRegex(ValueError,'채널 순서'):self.run_export()
        self.assertEqual(manifest.read_bytes(),before)

    def test_vcf_change_and_unknown_schema_fail_closed(self):
        import json
        from integrate_wgs_rna.pileup_features.export import validate_contig
        self.run_export()
        stamp=self.vcf.stat();os.utime(self.vcf,ns=(stamp.st_atime_ns,stamp.st_mtime_ns+1000000000))
        with self.assertRaisesRegex(ValueError,'지문'):self.run_export()
        path=self.directory/'contig_metadata.json';data=json.loads(path.read_text())
        data['pileup']['provenance']['schema']='unknown';path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,'schema'):validate_contig(self.directory)

    def test_independent_source_rejects_self_consistent_wrong_hp_content(self):
        import json,zarr
        from integrate_wgs_rna.pileup_features.export import validate_contig
        self.run_export()
        group=zarr.open_group(str(self.directory/'pileup.zarr'),mode='a')
        a=group['HAP1/features'];a[0,CI['hp1_read_A_ratio']]=.25
        path=self.directory/'contig_metadata.json';data=json.loads(path.read_text())
        data['pileup']['content_sha256']['HAP1']['features']=hashlib.sha256(a[:].tobytes()).hexdigest()
        path.write_text(json.dumps(data))
        validate_contig(self.directory,source=False)
        with self.assertRaisesRegex(ValueError,'독립 reference source'):validate_contig(self.directory,source=True)


if __name__ == '__main__':unittest.main()
