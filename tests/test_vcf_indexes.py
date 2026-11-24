"""누락index동시생성·원본보존·기존index보존·실패정리를검증함."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib,shutil,tempfile,unittest
import pysam
from handoff.vcf_indexes import ensure_index

class IndexTests(unittest.TestCase):
    def fixture(self,raw):
        p=Path(raw)/'a.vcf';p.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=40>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t2\t.\tC\tT\t.\tPASS\t.\n');q=Path(str(p)+'.gz');pysam.tabix_compress(str(p),str(q),force=True);return q
    def test_동시생성1회와원본불변(self):
        with tempfile.TemporaryDirectory() as raw:
            p=self.fixture(raw);before=p.read_bytes()
            with ThreadPoolExecutor(2) as pool:results=list(pool.map(lambda _:ensure_index(p,shutil.which('bcftools'),2),range(2)))
            self.assertCountEqual(results,['created','existing']);self.assertEqual(p.read_bytes(),before)
            with pysam.VariantFile(str(p)) as reader:self.assertEqual(len(list(reader.fetch('chr1'))),1)
    def test_기존인덱스보존(self):
        with tempfile.TemporaryDirectory() as raw:
            p=self.fixture(raw);pysam.tabix_index(str(p),preset='vcf');index=Path(str(p)+'.tbi');before=index.read_bytes()
            self.assertEqual(ensure_index(p,'실행되면실패'), 'existing');self.assertEqual(index.read_bytes(),before);self.assertFalse(Path(str(p)+'.csi').exists())
    def test_원본누락은생성하지않음(self):
        with tempfile.TemporaryDirectory() as raw:
            with self.assertRaises(FileNotFoundError):ensure_index(Path(raw)/'none.vcf.gz')
    def test_생성실패는부분인덱스남기지않음(self):
        import subprocess
        with tempfile.TemporaryDirectory() as raw:
            p=self.fixture(raw)
            with self.assertRaises(subprocess.CalledProcessError):ensure_index(p,'/bin/false')
            self.assertFalse(Path(str(p)+'.csi').exists());self.assertFalse(list(Path(raw).glob('.*.csi.*')))
if __name__=='__main__':unittest.main()
