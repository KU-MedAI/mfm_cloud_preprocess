"""GPU 통합 중복표시의 정책·cache·후속 단계 계약을 확인함."""
import importlib.util,json,os,shutil,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import pysam

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from integrate_wgs_rna import gpu_rna
from rna_cache import cache_valid
import test_rna_cache_contract as cache
from test_rna_cache_contract import write_bam,module,invoke
METRICS="READ_PAIRS_EXAMINED\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tUNPAIRED_READS_EXAMINED\tUNPAIRED_READ_DUPLICATES\tUNMAPPED_READS\tPERCENT_DUPLICATION\n1\t0\t0\t0\t0\t0\t0\n"


class IntegratedMarkdupTest(unittest.TestCase):
    def test_통합선택은RNA와BAM모두GPU일때만적용함(self):
        for rna,bam,mode,expected in [('cpu','parabricks','integrated',False),('parabricks','cpu','integrated',False),
                                      ('parabricks','parabricks','integrated',True),('parabricks','parabricks','separate',False)]:
            with patch.dict(os.environ,{'MFM_RNA_BACKEND':rna,'MFM_BAM_BACKEND':bam,'RNA_GPU_MARKDUP_MODE':mode}):
                self.assertEqual(gpu_rna.integrated_markdup_enabled(),expected)
        with patch.dict(os.environ,{'RNA_GPU_MARKDUP_MODE':'typo'}),self.assertRaises(ValueError):gpu_rna.integrated_markdup_enabled()

    def test_통합명령은광학100_coordinate정책과metrics를요청함(self):
        import rna_profile
        values=gpu_rna.parse(['--genomeDir','/index','--readFilesIn','/r1','/r2',*rna_profile.STAR_ALIGN_ARGS])
        with patch.dict(os.environ,{'MFM_GPU_RNA_POLICY':'experimental'}):
            command=gpu_rna.align_arguments(values,'/ref',Path('/out/a.bam'),4,74,integrated_markdup=True)
        self.assertNotIn('--no-markdups',command)
        self.assertNotIn('--markdups-assume-sortorder-queryname',command)
        self.assertEqual(command[command.index('--optical-duplicate-pixel-distance')+1],'100')
        self.assertEqual(command[command.index('--out-duplicate-metrics')+1],'/out/MarkDuplicates.metrics.txt')

    def test_중복표시단계는재계산없이복사_index_독립QC하고원본cache를보존함(self):
        stage=module('15')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ad=root/'aln';ad.mkdir();bam=ad/'star.bam';write_bam(bam)
            metrics=ad/'MarkDuplicates.metrics.txt'
            metrics.write_text(METRICS)
            record=gpu_rna.integrated_markdup_record(ad,bam)
            al=root/'align.json';al.write_text(json.dumps({'align_dir':str(ad),'star_bam':str(bam),'integrated_markdup':record}))
            calls=[]
            def run(cmd,**kwargs):
                calls.append(cmd)
                if 'index' not in cmd:raise AssertionError(cmd)
                pysam.index(cmd[-1])
            args=['--sample','S','--hap','1','--align-json',al,'--out-dir',root/'final',
                  '--work-dir',root,'--log-dir',root,'--picard',sys.executable,'--samtools',shutil.which('samtools')]
            with patch.dict(os.environ,{'MFM_BAM_BACKEND':'parabricks','MFM_CPU_BUDGET':'2'}),\
                 patch.object(stage,'run',run),patch.object(stage.parabricks_bam,'check_markdup_validation',wraps=stage.parabricks_bam.check_markdup_validation) as gate,\
                 patch.object(stage.parabricks_bam,'mark_duplicates',side_effect=AssertionError('GPU 재실행 금지')):
                invoke(stage,args);invoke(stage,args)
            self.assertEqual(len(calls),1)
            self.assertEqual(gate.call_count,2)
            self.assertTrue(cache_valid(record))
            result=json.loads((root/'hap1.markdup.json').read_text())
            self.assertTrue(result['backend_execution']['integrated_with_alignment'])
            self.assertTrue(result['resource_execution']['cache_reused'])
            self.assertEqual(result['primary_records'],2)
            self.assertNotEqual(bam.stat().st_ino,Path(result['marked_bam']).stat().st_ino)
            self.assertEqual(bam.read_bytes(),Path(result['final_bam']).read_bytes())
            metrics.write_text('changed')
            with self.assertRaisesRegex(ValueError,'계약/입력 변경'):gpu_rna.validate_integrated_markdup(record,bam)


class IntegratedCacheTest(cache.CacheContract):
    def fake_run(self,cmd,**kwargs):
        super().fake_run(cmd,**kwargs)
        if 'alignReads' in cmd:
            d=Path(cmd[cmd.index('--outFileNamePrefix')+1]);record=None
            if gpu_rna.integrated_markdup_enabled():
                (d/'MarkDuplicates.metrics.txt').write_text(METRICS)
                record=gpu_rna.integrated_markdup_record(d,d/'Aligned.sortedByCoord.out.bam')
            (d/'gpu_backend.json').write_text(json.dumps({'integrated_markdup':record}))

    def test_정렬cache는metrics손상과통합분리모드변경을감지함(self):
        self.build_index()
        with patch.dict(os.environ,{'MFM_RNA_BACKEND':'parabricks','MFM_BAM_BACKEND':'parabricks','RNA_GPU_MARKDUP_MODE':'integrated'}):
            first=self.build_align();second=self.build_align()
            self.assertEqual(first['align_dir'],second['align_dir'])
            self.assertTrue(first['integrated_markdup'])
            Path(first['integrated_markdup']['metrics']).write_text('changed')
            third=self.build_align();self.assertNotEqual(first['align_dir'],third['align_dir'])
            with patch.dict(os.environ,{'RNA_GPU_MARKDUP_MODE':'separate'}):
                fourth=self.build_align();self.assertIsNone(fourth['integrated_markdup'])
                self.assertNotEqual(third['align_dir'],fourth['align_dir'])


if __name__=='__main__':unittest.main()
