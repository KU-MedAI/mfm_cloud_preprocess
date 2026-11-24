"""실제BGZF검증증거의rename공유·재사용·변경검출·실패보존을검증함."""
from pathlib import Path
from unittest.mock import patch
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os,shutil,tempfile,time,unittest
import pysam
from handoff import integrity_cache as cache
from handoff import alignment_integrity as core
from handoff.downloads import fetch

class CacheTests(unittest.TestCase):
    def bam(self,path):
        with pysam.AlignmentFile(str(path),'wb',header={'HD':{'VN':'1.6'},'SQ':[{'SN':'chr1','LN':100}]}) as out:
            r=pysam.AlignedSegment(out.header);r.query_name='r';r.query_sequence='ACGT';r.reference_id=0;r.reference_start=0;r.cigarstring='4M';r.flag=0;out.write(r)
    def test_동일파일동시worker도한번만전수검사(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);path=root/'a.bam';self.bam(path);clock=time.time_ns()+3_000_000_000
            with patch.object(core,'scan',wraps=core.scan) as scan,patch.object(time,'time_ns',return_value=clock):
                with ThreadPoolExecutor(2) as pool:r=list(pool.map(lambda i:cache.validate_shared(path,root/f'worker{i}.json',root),range(2)))
                self.assertEqual(scan.call_count,1);self.assertTrue(all(x['status']=='pass' for x in r))
    def test_fetch의통제된rename검증증거를worker가재사용(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);source=root/'source.bam';self.bam(source);data=root/'data';data.mkdir();target=data/'a.bam';resource={'path':'a.bam','url':source.as_uri(),'kind':'wgs_alignment','size_bytes':source.stat().st_size}
            with patch.object(core,'scan',wraps=core.scan) as scan:
                fetch(resource,data,root/'state');self.assertEqual(scan.call_count,1)
                with patch.object(time,'time_ns',return_value=time.time_ns()+3_000_000_000):result=cache.validate_shared(target,root/'worker.json',data)
                self.assertEqual(scan.call_count,1);self.assertIn('relocated_from',result)
    def test_복사된다른inode는rename증거로게시불가(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);partial=root/'a.bam.part';self.bam(partial);payload=core.validate(partial,root/'proof.json');target=root/'a.bam';shutil.copyfile(partial,target)
            with self.assertRaisesRegex(ValueError,'다름'):cache.publish_relocated(payload,target,root/'final.json')
    def test_검증후손상은다시검출하고worker실패기록보존(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);path=root/'a.bam';self.bam(path)
            with patch.object(time,'time_ns',return_value=time.time_ns()+3_000_000_000):cache.validate_shared(path,root/'worker.json',root)
            rawdata=bytearray(path.read_bytes());block_size=int.from_bytes(rawdata[16:18],'little')+1;rawdata[block_size-8]^=3;path.write_bytes(rawdata)
            with self.assertRaises(core.AlignmentCorruption):cache.validate_shared(path,root/'worker.json',root)
            self.assertEqual(json.loads((root/'worker.json').read_text())['status'],'failed')
    def test_깨진공유JSON은전수검증으로복구(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);path=root/'a.bam';self.bam(path);state=root/'.download_records';state.mkdir();key=hashlib.sha256(str(path.resolve()).encode()).hexdigest();(state/(key+'.integrity.json')).write_text('{broken')
            self.assertEqual(cache.validate_shared(path,root/'worker.json',root)['status'],'pass')
    def test_미래mtime은대기없이전체검사유지(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);path=root/'a.bam';self.bam(path);future=time.time()+3600;os.utime(path,(future,future))
            with patch.object(core,'scan',wraps=core.scan) as scan,patch.object(time,'sleep',side_effect=AssertionError('미래mtime대기금지')):
                cache.validate_shared(path,root/'worker.json',root);cache.validate_shared(path,root/'worker.json',root);self.assertEqual(scan.call_count,2)
    def test_동일stat시간값의신규파일도CRC손상검출(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);path=root/'a.bam';self.bam(path)
            with patch.object(time,'time_ns',return_value=time.time_ns()+3_000_000_000):cache.validate_shared(path,root/'worker.json',root)
            prior=core.fingerprint(path);data=bytearray(path.read_bytes());size=int.from_bytes(data[16:18],'little')+1;data[size-8]^=1;path.write_bytes(data)
            with patch.object(core,'fingerprint',return_value=prior),patch.object(time,'time_ns',return_value=max(prior['mtime_ns'],prior['ctime_ns'])+1_000_000_000),patch.object(time,'sleep',side_effect=AssertionError('fresh검증을대기로대체하면안됨')):
                with self.assertRaises(core.AlignmentCorruption):cache.validate_shared(path,root/'worker.json',root)
if __name__=='__main__':unittest.main()
