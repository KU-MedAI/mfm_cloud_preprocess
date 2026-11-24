"""완료HAP보존·1회HAP2경계증설·잘못된요청차단을검증함."""
from pathlib import Path
from unittest.mock import patch
import json,os,tempfile,unittest,shutil
import pysam
from integrate_wgs_rna.realign.boundary_upgrade import prepare_boundary,finish_boundary,process_start
from integrate_wgs_rna.memory_policy import GIB

class BoundaryTests(unittest.TestCase):
    def fixture(self,raw):
        pool=Path(raw)/'pool';(pool/'leases').mkdir(parents=True);(pool/'waiters').mkdir()
        root=Path(raw)/'S';log=root/'logs/personal_realign';log.mkdir(parents=True)
        prior=root/'realigned/all_to_both/HAP1';prior.mkdir(parents=True)
        output=prior/'S.personal.HAP1.bam'
        with pysam.AlignmentFile(str(output),'wb',header={'HD':{'VN':'1.6'},'SQ':[{'SN':'chr1','LN':40}]}):pass
        pysam.index(str(output))
        metadata=prior/'S.realign.json';metadata.write_text(json.dumps({'sample':'S','haplotype':'HAP1','hap_read_mode':'all_to_both','format':'bam','qc_contract':{'status':'pass'}}))
        hap=root/'individual_fasta/HAP2';hap.mkdir(parents=True)
        for n in range(1,23):
            fasta=hap/f'chr{n}.fa';fasta.write_text(f'>chr{n}\n'+('ACGT'*10)+'\n');pysam.faidx(str(fasta))
        script=log/'S.HAP2.realign.sh';script.write_text('#!/bin/bash\ncat /dev/null | "$BWA_MEM2" mem -p -Y -C -t 2 index - | "$SAMTOOLS" sort -@ 0 -m 2G -o out -\n')
        lease=pool/'leases/owner.json';lease.write_text(json.dumps({'sample':'S','cpus':5,'total_cpus':16,'memory_bytes':80*GIB,'manager_pid':os.getpid(),'child_pid':None,'realign_phase':'align_restore_sort_qc','effective_alignment':{'hap_jobs':1}}))
        (pool/'leases/other.json').write_text(json.dumps({'cpus':4,'memory_bytes':4*GIB,'manager_pid':os.getpid(),'child_pid':None}))
        requests=pool/'hap_boundary_requests.json';requests.write_text(json.dumps({'samples':{'S':{'status':'armed','manager_pid':os.getpid(),'manager_start':process_start(os.getpid()),'sample_root':str(root),'hap_read_mode':'all_to_both'}}}))
        env={'MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_RAM_LIMIT_BYTES':str(500*GIB),'MFM_WGS_ALIGN_REQUESTED':'auto','MFM_WGS_SORT_REQUESTED':'auto','WGS_AUTO_THREADS':'1','SAMTOOLS_BIN':shutil.which('samtools')}
        return pool,root,script,log/'S.HAP2.resource.json',lease,requests,metadata,env

    def test_HAP1실제BAM검증후HAP2만증설하고1회적용(self):
        with tempfile.TemporaryDirectory() as raw:
            pool,root,script,report,lease,requests,metadata,env=self.fixture(raw);before=script.read_text();prior=metadata.read_bytes()
            with patch.dict(os.environ,env,clear=True):
                command,info=prepare_boundary(['bash',str(script)],report)
                self.assertGreater(info['threads'],2);self.assertEqual(command,['bash',str(script)])
                self.assertEqual(script.with_suffix('.before_boundary.sh').read_text(),before)
                self.assertEqual(metadata.read_bytes(),prior)
                self.assertEqual(json.loads(requests.read_text())['samples']['S']['status'],'applied')
                self.assertLessEqual(json.loads(lease.read_text())['cpus']+4,16)
                self.assertIsNone(prepare_boundary(command,report)[1])
                finish_boundary(info)
                row=json.loads(lease.read_text());self.assertEqual(row['realign_phase'],'output_qc');self.assertGreaterEqual(row['cpus'],5);self.assertLessEqual(row['cpus']+4,16)

    def test_PID시작시간다르면무수정(self):
        with tempfile.TemporaryDirectory() as raw:
            pool,root,script,report,lease,requests,metadata,env=self.fixture(raw);r=json.loads(requests.read_text());r['samples']['S']['manager_start']+=1;requests.write_text(json.dumps(r));before=script.read_bytes()
            with patch.dict(os.environ,env,clear=True):self.assertIsNone(prepare_boundary(['bash',str(script)],report)[1])
            self.assertEqual(script.read_bytes(),before)

    def test_경계sort버퍼와실제예약일치(self):
        with tempfile.TemporaryDirectory() as raw:
            pool,root,script,report,lease,requests,metadata,env=self.fixture(raw)
            (pool/'policy.json').write_text(json.dumps({'realign_sort_memory_bytes':16*GIB}))
            with patch.dict(os.environ,env,clear=True):
                _,info=prepare_boundary(['bash',str(script)],report)
            self.assertEqual(int(info['sort_mem']),16*GIB)
            self.assertIn('-m '+str(16*GIB),script.read_text())
            model=json.loads(lease.read_text())['effective_alignment']['memory_model']
            self.assertEqual(model['sort_buffers_bytes'],16*GIB*max(1,info['sort_threads']))

    def test_숫자thread명시한작업무수정(self):
        with tempfile.TemporaryDirectory() as raw:
            *_,env=self.fixture(raw);env['MFM_WGS_ALIGN_REQUESTED']='2';script=Path(raw)/'S/logs/personal_realign/S.HAP2.realign.sh'
            with patch.dict(os.environ,env,clear=True):self.assertIsNone(prepare_boundary(['bash',str(script)],script.with_name('S.HAP2.resource.json'))[1])

    def test_HAP1QC실패이면실행전차단(self):
        with tempfile.TemporaryDirectory() as raw:
            pool,root,script,report,lease,requests,metadata,env=self.fixture(raw);r=json.loads(metadata.read_text());r['qc_contract']['status']='fail';metadata.write_text(json.dumps(r));before=lease.read_bytes()
            with patch.dict(os.environ,env,clear=True):
                with self.assertRaisesRegex(ValueError,'HAP1'):prepare_boundary(['bash',str(script)],report)
            self.assertEqual(lease.read_bytes(),before)

    def test_요청없는일반실행은기존동작(self):
        with patch.dict(os.environ,{},clear=True):self.assertEqual(prepare_boundary(['echo','ok'],Path('/tmp/resource.json')),(['echo','ok'],None))

if __name__=='__main__':unittest.main()
