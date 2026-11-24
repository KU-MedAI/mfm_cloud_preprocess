"""GPU IO 선택 기능의 CPU 기본값·바이트/CRC/순서 보존 및 실패 차단을 검사함."""
from pathlib import Path
from unittest.mock import patch
import gzip,hashlib,io,json,os,random,struct,tempfile,unittest,zlib
import pysam
from handoff import acceleration,io_env,io_integrity,alignment_integrity
from integrate_wgs_rna import gpu_bgzf,gpu_nvcomp
from integrate_wgs_rna.gpu_io_stream import prefetch

class CPUCodec:
    def compress(self,chunks):
        result=[]
        for raw in chunks:
            c=zlib.compressobj(wbits=-15);result.append((c.compress(raw)+c.flush(),zlib.crc32(raw)))
        return result

class FrontendTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.bam=self.root/'source.bam'
        with pysam.AlignmentFile(str(self.bam),'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000000}]}) as f:
            for i in range(10):
                r=pysam.AlignedSegment();r.query_name='r'+str(i);r.query_sequence='ACGT'*20000;r.query_qualities=pysam.qualitystring_to_array('I'*80000);r.flag=1024 if i%2 else 0;r.reference_id=0;r.reference_start=i;r.mapping_quality=60;r.cigarstring='80000M';r.set_tag('HP',1);r.set_tag('PS',5);f.write(r)
    def test_CPU기본환경_설정불변(self):
        original={'pipeline_env':{'WINDOW_STORAGE_POLICY':'annotated_acgtn_v1'}}
        with patch.dict(os.environ,{},clear=True):self.assertEqual(acceleration.apply(original),original)
    def test_IO명시선택은_독립설정임(self):
        c=acceleration.apply({'pipeline_env':{}},{'io_backend':'nvcomp'})
        self.assertEqual(c['pipeline_env']['MFM_IO_BACKEND'],'nvcomp')
        self.assertNotIn('MFM_WGS_BACKEND',c['pipeline_env'])
    def test_CPU명시가_낡은GPU설정을해제(self):
        with patch.dict(os.environ,{'MFM_IO_BACKEND':'nvcomp'}):
            io_env.apply({'acceleration':{'io_backend':'cpu'}});self.assertEqual(os.environ['MFM_IO_BACKEND'],'cpu')
    def test_다운로드에도_허용GPU와_IO설정전달(self):
        c={'pipeline_env':{'MFM_IO_BACKEND':'nvcomp','MFM_GPU_DEVICES':'0,1,2','CUDA_VISIBLE_DEVICES':'0,1,2'}}
        with patch.dict(os.environ,{},clear=True):
            io_env.apply(c);self.assertEqual(os.environ['MFM_GPU_DEVICES'],'0,1,2');self.assertEqual(os.environ['MFM_IO_BACKEND'],'nvcomp')
    def test_BGZF전체해독바이트_SHA_순서보존(self):
        output=io.BytesIO()
        with self.bam.open('rb') as f:report=gpu_bgzf.compress_stream(f,output,CPUCodec(),batch_chunks=2,verify_workers=2)
        raw=gzip.decompress(self.bam.read_bytes());self.assertEqual(gzip.decompress(output.getvalue()),raw)
        self.assertEqual(report['raw_sha256'],hashlib.sha256(raw).hexdigest())
        p=self.root/'output.bam';p.write_bytes(output.getvalue());alignment_integrity.scan(p);pysam.index(str(p))
        with pysam.AlignmentFile(str(self.bam),'rb') as a,pysam.AlignmentFile(str(p),'rb') as b:self.assertEqual([r.to_string() for r in a],[r.to_string() for r in b])
    def test_CRC불일치_무조건차단(self):
        class Bad(CPUCodec):
            def compress(self,chunks):return [(x,c^1) for x,c in super().compress(chunks)]
        with self.bam.open('rb') as f,self.assertRaises(ValueError):gpu_bgzf.compress_stream(f,io.BytesIO(),Bad(),batch_chunks=2)
    def test_GPU압축출력손상_무조건차단(self):
        class Bad(CPUCodec):
            def compress(self,chunks):return [(b'bad',c) for x,c in super().compress(chunks)]
        with self.bam.open('rb') as f,self.assertRaises((ValueError,zlib.error)):gpu_bgzf.compress_stream(f,io.BytesIO(),Bad(),batch_chunks=2)
    def test_GPU출력배치누락_차단(self):
        class Bad(CPUCodec):
            def compress(self,chunks):return super().compress(chunks)[:-1]
        with self.bam.open('rb') as f,self.assertRaises(ValueError):gpu_bgzf.compress_stream(f,io.BytesIO(),Bad(),batch_chunks=2)
    def test_EOF누락_차단(self):
        with self.assertRaises(ValueError):gpu_bgzf.compress_stream(io.BytesIO(self.bam.read_bytes()[:-28]),io.BytesIO(),CPUCodec())
    def test_원본BAM검사는_CPU검출기유지(self):
        with patch.dict(os.environ,{'MFM_IO_BACKEND':'nvcomp'}),patch.object(alignment_integrity,'validate',return_value={'status':'pass'}) as check:
            self.assertEqual(io_integrity.validate(self.bam,self.root/'physical.json')['status'],'pass');check.assert_called_once()
    def test_GPU검사_잘못된형식의_ValueError보존(self):
        from types import SimpleNamespace
        path=self.root/'unsupported.cram';path.write_bytes(b'CRAM'+bytes([4,0])+bytes(30))
        def child(command,**kwargs):
            report=Path(command[command.index('--report')+1]);report.write_text(json.dumps({'status':'failed','error':'unsupported version','error_type':'ValueError','corruption':False}))
            return SimpleNamespace(returncode=1,stderr='')
        with patch.dict(os.environ,{'MFM_IO_BACKEND':'nvcomp','MFM_GPU_CRAM_MIN_BYTES':'0'},clear=True),patch.object(io_integrity.subprocess,'run',side_effect=child):
            with self.assertRaisesRegex(ValueError,'unsupported version'):io_integrity.validate(path,self.root/'physical.json')
    def test_bounded_prefetch_예외전달(self):
        def data():
            yield 1;yield 2;raise ValueError('fixture error')
        iterator=prefetch(data(),buffers=1);self.assertEqual(next(iterator),1);self.assertEqual(next(iterator),2)
        with self.assertRaises(ValueError):next(iterator)
    def test_크기상한_차단(self):
        with self.assertRaises(ValueError):gpu_bgzf.compress_stream(io.BytesIO(),io.BytesIO(),CPUCodec(),batch_chunks=2000)
    def test_미지원_IO설정차단(self):
        with self.assertRaises(ValueError):acceleration.apply({}, {'io_backend':'bad'})

if __name__=='__main__':unittest.main()
