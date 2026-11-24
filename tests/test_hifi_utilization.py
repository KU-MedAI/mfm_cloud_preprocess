"""실제HiFi재정렬·입력BAM/CRAM집계·동시HAP예약·단계투입을검증함."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from dataclasses import replace
import json,os,random,shutil,tempfile,unittest
import pysam
from integrate_wgs_rna.realign.config import PersonalRealignConfig
from integrate_wgs_rna.realign.qc import _count_input_reads_one_pass
from integrate_wgs_rna.realign.phases import RealignPhases,_atomic
from integrate_wgs_rna.realign.hap_resources import HaplotypeResources
from integrate_wgs_rna.realign.align_resources import alignment_plan
from integrate_wgs_rna.personal_realign import run_personal_realign_generation
from integrate_wgs_rna.memory_policy import GIB
from handoff.stage_resources import longread_memory_plan,starting_memory
from handoff.queue_resources import active_ram,full_ram

class IndexedCounts(unittest.TestCase):
    def test_BAM_CRAM_배치미배치unmapped_HP이상값을기존과비교(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);ref=root/'ref.fa';ref.write_text('>chr1\n'+'A'*200+'\n>chr2\n'+'A'*200+'\n');pysam.faidx(str(ref));bam=root/'a.bam'
            with pysam.AlignmentFile(str(bam),'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':200},{'SN':'chr2','LN':200}]}) as out:
                for tid in [0,1]:
                    for i,hp in enumerate([None,0,1,2,9,'bad']):
                        r=pysam.AlignedSegment(out.header);r.query_name=f'r{tid}_{i}';r.query_sequence='AAAA';r.query_qualities=pysam.qualitystring_to_array('IIII');r.reference_id=tid;r.reference_start=i*10;r.flag=0;r.cigarstring='4M';r.mapping_quality=60
                        if hp is not None:r.set_tag('HP',hp)
                        out.write(r)
                    r=pysam.AlignedSegment(out.header);r.query_name=f'placed{tid}';r.query_sequence='AAAA';r.reference_id=tid;r.reference_start=100;r.flag=4;out.write(r)
                for i in range(4):
                    r=pysam.AlignedSegment(out.header);r.query_name=f'unplaced{i}';r.query_sequence='AAAA';r.reference_id=-1;r.reference_start=-1;r.flag=4;out.write(r)
            pysam.index(str(bam));cram=root/'a.cram'
            with pysam.AlignmentFile(str(bam)) as src,pysam.AlignmentFile(str(cram),'wc',header=src.header,reference_filename=str(ref),format_options=[b'seqs_per_slice=1',b'slices_per_container=1']) as dst:
                for r in src:dst.write(r)
            pysam.index(str(cram));cfg=SimpleNamespace(reference_fasta=ref,sort_threads=0,hap_jobs=1,jobs=1,samtools_bin=shutil.which('samtools'))
            for path in [bam,cram]:
                values=[]
                for backend in ['serial','indexed']:
                    with patch.dict(os.environ,{'MFM_REALIGN_INPUT_BACKEND':backend,'MFM_REALIGN_INPUT_CHUNK_BP':'15','MFM_CPU_BUDGET':'4','MFM_SAMPLE_CPU_BUDGET':'4','MFM_RAM_LIMIT_BYTES':str(200*GIB),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'},clear=True):values.append(_count_input_reads_one_pass(cfg,path,['chr1'],invalid_hp_as_untagged=False))
                self.assertEqual(values[0],values[1]);self.assertEqual(values[1],{'hp1':1,'hp2':1,'untagged':2,'all_primary_mapped':6,'unmapped':6})

class ResourceTests(unittest.TestCase):
    def test_다른HAP정렬예약은반환하지않고QC만자원변경(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();(root/'waiters').mkdir();path=root/'leases/owner.json'
            _atomic(path,{'manager_pid':os.getpid(),'child_pid':None,'cpus':24,'total_cpus':32,'memory_bytes':80*GIB});_atomic(root/'leases/other.json',{'manager_pid':os.getpid(),'child_pid':None,'cpus':8,'memory_bytes':8*GIB})
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path),'MFM_CPU_BUDGET':'32','MFM_RAM_LIMIT_BYTES':str(500*GIB)},clear=True):
                group=HaplotypeResources(RealignPhases(),None,[1,2]);a=group.controller(1);b=group.controller(2)
                self.assertEqual(a.reserve('index',max_units=5,base_memory=8*GIB),5)
                row=json.loads(path.read_text());self.assertEqual(row['cpus'],17);self.assertEqual(group.states[2]['cpus'],12);self.assertEqual(group.states[2]['memory_bytes'],40*GIB)
                group.release(1);self.assertEqual(json.loads(path.read_text())['cpus'],12)
                b.reserve('QC',max_units=20,memory_per_unit=4*GIB);self.assertLessEqual(json.loads(path.read_text())['cpus']+8,32)
                self.assertEqual(b.resource_limits['cpus'],20);group.release(2);self.assertEqual(json.loads(path.read_text())['cpus'],1)

    def test_minimap보조IO와unmappedawk도CPU예약에포함(self):
        c=SimpleNamespace(platform='hifi',threads=9,sort_threads=3,hap_jobs=2,sort_mem='2G')
        p=alignment_plan(c,{1:3000000000,2:3000000000},25,900000000000,{'WGS_AUTO_THREADS':'1','MFM_WGS_ALIGN_REQUESTED':'auto','MFM_WGS_SORT_REQUESTED':'auto','MFM_WGS_HAP_REQUESTED':'auto'})
        self.assertLessEqual(p['reserved_cpus'],25);self.assertEqual(p['reserved_cpus'],p['hap_jobs']*(p['threads']+p['sort_threads']+5+p['input_decode_threads']))

    def test_새후보는64GiB로시작하고미래피크검사는보존(self):
        row={'sample_id':'S','platform':'hifi','read_type':'long_read','vcf_source':'public','phase_source':'provided','wgs_input_format':'bam','data_group':'wgs'}
        self.assertEqual(starting_memory(row),64*GIB);self.assertGreater(longread_memory_plan(row)['admission_bytes'],starting_memory(row))
        with patch('handoff.queue_resources.live_leases',return_value={}):self.assertEqual(active_ram(Path('/tmp'),[],row),16*GIB+64*GIB)
        self.assertEqual(starting_memory(row,{'pipeline_env':{'MFM_STAGE_ADMISSION':'0'}}),longread_memory_plan(row)['admission_bytes'])
        self.assertEqual(starting_memory(row,{'sample_pipeline_env':{'S':{'WGS_START_AT':'10'}}}),4*GIB)
        self.assertGreater(starting_memory(row,{'pipeline_env':{'WGS_START_AT':'10','MFM_REALIGN_PHASE_RESOURCES':'0'}}),4*GIB)

class ActualHiFi(unittest.TestCase):
    def test_실제minimap_직렬병렬단계재배분후read_QC동일(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);rng=random.Random(71);sequence=''.join(rng.choice('ACGT') for _ in range(8000));ref=root/'ref.fa';ref.write_text('>chr1\n'+sequence+'\n');pysam.faidx(str(ref));source=root/'source.bam'
            with pysam.AlignmentFile(str(source),'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':len(sequence)}],'RG':[{'ID':'r','SM':'S','PL':'PACBIO'}]}) as out:
                for i in range(8):
                    r=pysam.AlignedSegment(out.header);r.query_name=f'r{i}';r.query_sequence=sequence[i*700:i*700+1000];r.query_qualities=pysam.qualitystring_to_array('I'*1000);r.reference_id=0;r.reference_start=i*700;r.flag=0;r.mapping_quality=60;r.cigarstring='1000M';r.set_tag('HP',i%3);r.set_tag('PS',1);r.set_tag('RG','r');out.write(r)
            pysam.index(str(source));hapdirs={}
            for hap in [1,2]:
                folder=root/f'inputHAP{hap}';folder.mkdir();f=folder/'chr1.fa';f.write_text('>chr1\n'+sequence+'\n');pysam.faidx(str(f));(folder/'chr1.chain').write_text(f'chain 0 chr1 8000 + 0 8000 chr1 8000 + 0 8000 1\n8000\n');hapdirs[hap]=folder
            outputs=[]
            for label,phase in [('old','0'),('new','1'),('pair16','1'),('pair32','1')]:
                cfg=PersonalRealignConfig(cohort='T',output_root=root/label,sample_list=None,samples=['S'],chroms=['chr1'],sequence_variant_policy='snv_indel',dataset_role='pre_training',platform='hifi',realigned_subdir='realigned/all_to_both',alignment=source,hap1_fasta_dir=hapdirs[1],hap2_fasta_dir=hapdirs[2],reference_fasta=ref,hap_read_mode='all_to_both',threads=2,sort_threads=1,sort_mem='32M',hap_jobs=2,include_unmapped=False,inspect_reads=5)
                env={'PATH':os.environ['PATH'],'MFM_REALIGN_PHASE_RESOURCES':phase,'MFM_CPU_BUDGET':'24','MFM_SAMPLE_CPU_BUDGET':'24','MFM_RAM_LIMIT_BYTES':str(500*GIB),'MFM_REALIGN_INPUT_BACKEND':'indexed','MFM_WGS_ALIGN_REQUESTED':'auto','MFM_WGS_SORT_REQUESTED':'auto','MFM_WGS_HAP_REQUESTED':'auto','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
                env['PYTHONPATH']=os.environ.get('PYTHONPATH',str(Path(__file__).resolve().parents[1]/'src'))
                if label in {'pair16','pair32'}:
                    cfg.sort_mem='2G'
                    cpus='16' if label=='pair16' else '32'
                    env.update(MFM_CPU_BUDGET=cpus, MFM_SAMPLE_CPU_BUDGET=cpus,
                               MFM_WGS_SORT_MEM_REQUESTED='auto')
                if phase=='1':
                    pool=root/f'pool_{label}';(pool/'leases').mkdir(parents=True);(pool/'waiters').mkdir();lease=pool/'leases/owner.json';_atomic(lease,{'manager_pid':os.getpid(),'child_pid':None,'cpus':4,'total_cpus':int(env['MFM_CPU_BUDGET']),'memory_bytes':4*GIB});env['MFM_CPU_POOL_LEASE_FILE']=str(lease)
                with patch.dict(os.environ,env,clear=True):run_personal_realign_generation(cfg)
                result=[]
                for hap in [1,2]:
                    output=root/label/f'S/realigned/all_to_both/HAP{hap}/S.personal.HAP{hap}.cram';target=output.parent/'target.fa'
                    with pysam.AlignmentFile(str(output),reference_filename=str(target)) as bam:result.append([r.to_string() for r in bam])
                    metadata=json.loads((output.parent/'S.realign.json').read_text());self.assertEqual(metadata['qc_contract']['status'],'pass')
                outputs.append(result)
                if phase=='1':
                    self.assertEqual(json.loads(lease.read_text())['cpus'],1)
                    timing=json.loads((root/label/'S/qc/S.stage11_timing.json').read_text())
                    if timing['effective_threads']['hap_jobs']>1:
                        self.assertEqual(set(timing['hap_resource_phases']),{'1','2'})
                    else:
                        self.assertEqual(set(timing['effective_threads_by_hap']),{'HAP1','HAP2'})
                    self.assertGreater(timing['effective_threads']['input_decode_per_hap'],0)
                    self.assertIn('view -@', (root/label/'S/logs/personal_realign/S.HAP1.realign.sh').read_text())
                    if label in {'pair16','pair32'}:
                        self.assertEqual(timing['effective_threads']['hap_jobs'], 1 if label=='pair16' else 2)
                        for hap in [1,2]:
                            script=(root/label/f'S/logs/personal_realign/S.HAP{hap}.realign.sh').read_text()
                            self.assertIn('-t 10', script)
                            self.assertIn('-m '+str(8*GIB), script)
            for output in outputs[1:]:
                self.assertEqual(outputs[0],output)
if __name__=='__main__':unittest.main()
