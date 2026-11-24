"""캐시 입력 동일성·좌표·용량 상한·프로세스 재연결을 검증함."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json,multiprocessing,os,pickle,sys,tempfile,unittest
from unittest.mock import patch
import numpy as np
import zarr
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from integrate_wgs_rna.window_reader import WindowReader
from integrate_wgs_rna.model_input import load_window
from integrate_wgs_rna.window_chunk_cache import DecodedChunkCache


def _worker_read(reader):
    value=reader.load('chr1','HAP1',0)
    return reader.cache_info()['pid'],value['dna_codes'].tolist()


class ReaderTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)/'S';r=self.root
        (r/'qc').mkdir(parents=True);(r/'qc/final_package_validation.json').write_text(json.dumps({'sample':'S','status':'pass','failed_checks':[]}))
        layout=r/'compact/layout';layout.mkdir(parents=True);(layout/'manifest.json').write_text(json.dumps({'status':'complete','validation_level':'full'}))
        for contig in ['chr1','chr2']:
            seq=zarr.open_group(str(layout/contig/'hap_sequence.zarr'),mode='w');seq.attrs['sample']='S'
            dep=zarr.open_group(str(layout/contig/'depth.zarr'),mode='w')
            for hap in ['HAP1','HAP2']:
                g=seq.create_group(hap)
                codes=np.array([0,4,1,2,3,5,5,5]+[0,1,2,3]*2,dtype=np.uint8)
                for k,v in {'bases_packed':(codes[::2]<<4)|codes[1::2],'offsets':np.array([0,4,8]),'lengths':np.array([8,8]),'actual_length':np.array([5,8]),'raw_start':np.array([10,20])}.items():g.array(k,v)
                g=dep.create_group(hap);g.array('offsets',np.array([0,8,16]));g.array('depth',np.arange(16,dtype=np.uint32))
        for hap in ['HAP1','HAP2']:
            g=zarr.open_group(str(r/'rna_personal/panel'/hap/'coverage_windows.zarr'),mode='w');g.attrs.update(sample='S',hap=hap)
            for k,v in {'contig':np.array(['chr2','chr1','chr1','chr2']),'personal_start':np.array([20,10,20,10]),'personal_end':np.array([28,15,28,15]),'offsets':np.array([0,8,13,21,26]),'depth':np.arange(26,dtype=np.uint32)}.items():g.array(k,v)

    def tearDown(self):self.temp.cleanup()

    def test_same_values_reordered_rna_lru_and_padding(self):
        reader=WindowReader(self.root,max_cached_contigs=1)
        with patch('integrate_wgs_rna.window_reader.zarr.open_group',wraps=zarr.open_group) as opened:
            for _ in range(2):
                actual=reader.load('chr1','HAP1',0)
                np.testing.assert_array_equal(actual['rna_depth'],[8,9,10,11,12,0,0,0])
                self.assertEqual(actual['wgs_depth'][1],1)
            self.assertEqual(opened.call_count,3)
        for c in ['chr1','chr2','chr1']:
            for h in ['HAP1','HAP2']:
                for i in [0,1]:
                    a=reader.load(c,h,i);b=load_window(self.root,c,h,i)
                    for key in a:
                        if isinstance(a[key],np.ndarray):np.testing.assert_array_equal(a[key],b[key],err_msg=key)
                        else:self.assertEqual(a[key],b[key])
            self.assertEqual(reader.cache_info()['cached_contigs'],1)
            self.assertLessEqual(reader.cache_info()['estimated_metadata_bytes'],64<<20)
        reader.close();self.assertEqual(reader.cache_info()['estimated_metadata_bytes'],0)

    def test_sequence_only_preserves_rna_and_masks(self):
        import shutil
        layout=self.root/'compact/layout'
        (layout/'manifest.json').write_text(json.dumps({'status':'complete','validation_level':'full','wgs_payload':'sequence_only'}))
        for contig in ['chr1','chr2']:shutil.rmtree(layout/contig/'depth.zarr')
        reader=WindowReader(self.root)
        value=reader.load('chr1','HAP1',0)
        self.assertIsNone(value['wgs_depth'])
        np.testing.assert_array_equal(value['rna_depth'],[8,9,10,11,12,0,0,0])
        np.testing.assert_array_equal(value['dna_codes'],[0,4,1,2,3,5,5,5])
        self.assertEqual(reader.window_count('chr1','HAP1'),2)
        reader.close()

    def test_qc_coordinate_mutation_duplicates_and_budget(self):
        reader=WindowReader(self.root);reader.load('chr1','HAP1',0)
        p=self.root/'compact/layout/chr1/hap_sequence.zarr/HAP1/raw_start/0';p.touch()
        s=p.stat();os.utime(p,ns=(s.st_atime_ns,s.st_mtime_ns+1000000000))
        with self.assertRaisesRegex(ValueError,'metadata 변경'):reader.load('chr1','HAP1',0)
        with self.assertRaisesRegex(ValueError,'예산'):WindowReader(self.root,max_metadata_bytes=32).load('chr1','HAP1',0)
        g=zarr.open_group(str(self.root/'rna_personal/panel/HAP1/coverage_windows.zarr'),mode='a')
        g['contig'][2]='chr1';g['personal_start'][2]=10;g['personal_end'][2]=15
        with self.assertRaisesRegex(ValueError,'중복'):WindowReader(self.root).load('chr1','HAP1',0)
        q=self.root/'qc/final_package_validation.json';q.write_text(json.dumps({'sample':'S','status':'fail'}))
        with self.assertRaisesRegex(ValueError,'QC/manifest'):reader.load('chr1','HAP1',0)

    def test_decoded_chunk_lru_budget_and_change_detection(self):
        path=self.root/'fixture.zarr';array=zarr.open_array(str(path),mode='w',shape=(256,),chunks=(64,),dtype='u4')
        array[:]=np.arange(256,dtype=np.uint32)
        cache=DecodedChunkCache(256)
        np.testing.assert_array_equal(cache.read(array,0,30),np.arange(30))
        np.testing.assert_array_equal(cache.read(array,2,6),np.arange(2,6));self.assertEqual(cache.hits,1)
        np.testing.assert_array_equal(cache.read(array,60,80),np.arange(60,80));self.assertEqual(cache.bytes,256)
        self.assertEqual(len(cache.entries),1)
        p=path/'1';st=p.stat();os.utime(p,ns=(st.st_atime_ns,st.st_mtime_ns+1000000000))
        with self.assertRaisesRegex(ValueError,'chunk 변경'):cache.read(array,70,80)
        zero=DecodedChunkCache(0);np.testing.assert_array_equal(zero.read(array,0,10),array[:10]);self.assertEqual(zero.bytes,0)

    def test_bounded_axis_audit_with_actual_reader(self):
        import csv
        from integrate_wgs_rna.input_audit import audit_sample
        d=self.root/'metadata';d.mkdir()
        with (d/'wgs_windows.tsv').open('w') as f:
            fields=['contig','haplotype','window_index','raw_start','actual_length','n_count','n_fraction','kept']
            writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t');writer.writeheader()
            for c in ['chr1','chr2']:
                for h in ['HAP1','HAP2']:
                    for i in [0,1]:writer.writerow(dict(zip(fields,[c,h,i,10 if i==0 else 20,5 if i==0 else 8,1 if i==0 else 0,0.2 if i==0 else 0,1])))
        result=audit_sample(self.root,Path(self.temp.name)/'audit.json',workers=2,max_per_axis=1)
        self.assertEqual(result['status'],'pass');self.assertEqual(result['totals']['windows'],4)
        self.assertTrue(all(a['windows']==1 for a in result['axes']))
        self.assertFalse((self.root/'audit.json').exists())

    def test_spawn_pickle_and_fork_reconnect(self):
        reader=WindowReader(self.root);reader.load('chr1','HAP1',0)
        fresh=pickle.loads(pickle.dumps(reader));self.assertEqual(fresh.cache_info()['cached_contigs'],0)
        for method in ['spawn','fork']:
            if method not in multiprocessing.get_all_start_methods():continue
            with ProcessPoolExecutor(max_workers=1,mp_context=multiprocessing.get_context(method)) as pool:
                pid,codes=pool.submit(_worker_read,reader).result()
            self.assertNotEqual(pid,os.getpid());self.assertEqual(codes,[0,4,1,2,3,5,5,5])


if __name__=='__main__':unittest.main()
