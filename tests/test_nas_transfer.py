"""실물 소형 파일·rsync로 수동 이관/검증/회수/재시작 계약을 검사함."""
from pathlib import Path
from types import SimpleNamespace
import copy
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src'),str(ROOT/'vendor/MFM-DNARNAFM')]
from handoff import cli,nas_transfer as transfer,archive_receipt,local_queue
from integrate_wgs_rna.pileup_features.schema import FEATURE_SET

spec=importlib.util.spec_from_file_location('nas_receiver',ROOT/'tools/pull_completed.py')
receiver=importlib.util.module_from_spec(spec);spec.loader.exec_module(receiver)


class NasTransferTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='nas transfer ');self.addCleanup(self.temp.cleanup)
        self.base=Path(self.temp.name);self.state=self.base/'state';self.output=self.base/'source cohort'
        self.nas=self.base/'NAS/final_pileup';self.nas.parent.mkdir(parents=True)
        (self.state/'samples').mkdir(parents=True)
        self.config=dict(state_root=str(self.state),output_root=str(self.output),data_root=str(self.base/'raw'),
            run_id='test',selection_fingerprint='a'*64,model_target_bytes=10**12,cleanup_shared_fasta=False,
            cpu_pool_limit=4,memory_limit_bytes=10**12,sample_cpus=2,local_jobs=2,download_jobs=1,
            pipeline_env={'MFM_SHORT_READ_OUTPUT_MODE':'pileup_feature','MFM_RNA_OUTPUT_MODE':'full',
                          'WINDOW_STORAGE_POLICY':'annotated_acgtn_v1','CLEAN_INTERMEDIATES_ON_SUCCESS':'1'})
        self.rows=[]
        for sample in ['S1','S2']:
            row=dict(sample_id=sample,data_group='wgs_rna_pair',project_id='P',cohort='P',
                     platform='illumina',read_type='short_read',rna_path='rna/'+sample+'.bam',
                     known_input_total_bytes='1000',estimated_model_bytes='8')
            self.rows.append(row);root=self.sample(sample)
            for folder in ['metadata','qc','compact/fake/chr1','rna_personal/HAP1/coverage.zarr','logs']:
                (root/folder).mkdir(parents=True)
            (root/'compact/fake/chr1/0').write_bytes(b'ACGT')
            (root/'rna_personal/HAP1/coverage.zarr/0').write_bytes(b'0001')
            transfer.save(root/'metadata/sample.json',dict(sample_id=sample,layout_name='fake',
                sequence_variant_policy='snv_indel',dataset_role='post_training',read_type='short_read',
                project_id='P',hap_read_mode='all_to_both'))
            transfer.save(root/'metadata/window_policy.json',{'storage_policy':'annotated_acgtn_v1'})
            transfer.save(root/'metadata/cleanup_completed.json',dict(schema='mfm_delivery_cleanup_v1',
                sample=sample,mode='rna',dry_run=False,status='complete',post_cleanup_validation={
                    'schema':'mfm_delivery_final_validation_v1','sample':sample,'mode':'rna','status':'pass','errors':[]}))
            transfer.save(root/'qc/final_package_validation.json',dict(schema='mfm_final_package_validation_v1',
                sample=sample,mode='rna',status='pass',failed_checks=[],wgs_payload='pileup_feature',
                rna_payload='full',pileup_feature_set=FEATURE_SET))
            transfer.save(self.state/'samples'/(sample+'.json'),dict(sample=sample,status='success',model_bytes=8))
        from handoff.planner import write_tsv
        write_tsv(self.state/'selected_samples.tsv',self.rows)
        transfer.save(self.state/'effective_config.json',self.config)

    def sample(self,name):return self.output/'wgs_rna_pair/P'/name

    def prepare(self,**options):
        return transfer.prepare(self.config,self.rows,str(self.nas),**options)

    def receipt(self,batch,index=0):
        item=batch['samples'][index]
        return {**{k:item[k] for k in ['sample','destination','inventory_sha256','file_count','total_bytes']},
            'batch_id':batch['batch_id'],'token':batch['token'],'verification':'sha256_all_files','receiver':'test NAS'}

    def archive(self):
        batch=self.prepare(batch_size=1)
        result=transfer.commit(self.config,self.rows,self.receipt(batch),True)
        return batch,result

    def test_default_waits_for_ten_without_hashing_or_deleting(self):
        with patch.object(transfer,'inventories') as scan:
            self.assertEqual(self.prepare(),{'status':'waiting','ready':2,'required':10})
            scan.assert_not_called()
        self.assertTrue((self.sample('S1')/'compact').is_dir())

    def test_partial_tail_and_same_batch_resume(self):
        batch=self.prepare(allow_partial=True)
        self.assertEqual(len(batch['samples']),2)
        self.assertEqual(self.prepare()['batch_id'],batch['batch_id'])
        with self.assertRaises(ValueError):
            transfer.prepare(self.config,self.rows,'/tmp/other/final_pileup')
        with self.assertRaises(ValueError):
            transfer.prepare({**self.config,'run_id':'different'},self.rows,str(self.nas),batch_id=batch['batch_id'])

    def test_final_directory_is_separate_and_existing_final_untouched(self):
        with self.assertRaisesRegex(ValueError,'final_pileup'):
            transfer.prepare(self.config,self.rows,str(self.nas.parent/'final'))
        sentinel=self.nas.parent/'final/short_read_sequence/keep';sentinel.parent.mkdir(parents=True);sentinel.write_text('old')
        self.archive();self.assertEqual(sentinel.read_text(),'old')

    def test_destination_is_classified_by_read_group_project_and_sample(self):
        for read_type,folder in [('short_read','short_read_pileup'),('long_read','long_read_pileup')]:
            for group,subfolder in [('wgs','wgs_only'),('wgs_rna_pair','pair')]:
                row={**self.rows[0],'read_type':read_type,'data_group':group}
                with self.subTest(read_type=read_type,group=group):
                    self.assertEqual(transfer.destination_relative(row),f'{folder}/{subfolder}/P/S1')
        batch=self.prepare(batch_size=1);item=batch['samples'][0]
        self.assertEqual(item['relative_destination'],'short_read_pileup/pair/P/S1')
        self.assertEqual(item['destination'],str(self.nas/'short_read_pileup/pair/P/S1'))

    def test_unknown_classification_and_nested_project_are_rejected(self):
        for change in [{'read_type':''},{'read_type':'unknown'},{'data_group':'other'},
                       {'project_id':'P/subfolder'},{'project_id':'..'},{'project_id':'.'}]:
            with self.subTest(change=change),self.assertRaises(ValueError):
                transfer.destination_relative({**self.rows[0],**change})

    def test_prepared_batch_keeps_its_recorded_destination_on_resume(self):
        batch=self.prepare(batch_size=1)
        item=batch['samples'][0];item['relative_destination']='old_run/wgs_rna_pair/P/S1'
        item['destination']=str(self.nas/item['relative_destination'])
        transfer.save(transfer.batch_path(self.config,batch['batch_id']),batch)
        self.assertEqual(self.prepare()['samples'][0]['destination'],item['destination'])

    def test_not_completed_locked_or_wrong_policy_samples_are_not_selected(self):
        import fcntl
        transfer.save(self.state/'samples/S2.json',dict(sample='S2',status='running',model_bytes=0))
        with (self.state/'samples/S1.lock').open('a') as handle:
            fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
            self.assertEqual(self.prepare()['ready'],0)
        transfer.save(self.sample('S1')/'metadata/window_policy.json',{'storage_policy':'different'})
        self.assertEqual(self.prepare()['ready'],0)

    def test_wrong_receipt_never_deletes_source(self):
        batch=self.prepare(batch_size=1)
        for key,value in [('token','bad'),('destination','/tmp/final_pileup/wrong'),
                          ('file_count',999),('inventory_sha256','0'*64),('verification','size_only')]:
            with self.subTest(key=key),self.assertRaises(ValueError):
                transfer.commit(self.config,self.rows,{**self.receipt(batch),key:value},True)
            self.assertTrue((self.sample('S1')/'compact/fake/chr1/0').exists())

    def test_source_content_change_with_same_size_and_mtime_is_rejected(self):
        batch=self.prepare(batch_size=1);path=self.sample('S1')/'compact/fake/chr1/0';stat=path.stat()
        path.write_bytes(b'TGCA');os.utime(path,ns=(stat.st_atime_ns,stat.st_mtime_ns))
        with self.assertRaisesRegex(ValueError,'source 파일 변경'):
            transfer.commit(self.config,self.rows,self.receipt(batch),True)
        self.assertTrue(path.exists())

    def test_new_files_and_inventory_changes_are_rejected(self):
        batch=self.prepare(batch_size=1);(self.sample('S1')/'compact/new').write_bytes(b'new')
        with self.assertRaises(ValueError):transfer.commit(self.config,self.rows,self.receipt(batch),True)
        Path(batch['samples'][0]['inventory_path']).write_text('{}\n')
        with self.assertRaisesRegex(ValueError,'목록 변경'):
            transfer.commit(self.config,self.rows,self.receipt(batch),True)

    def test_copy_only_keeps_payload_and_records_verified_location(self):
        batch=self.prepare(batch_size=1);result=transfer.commit(self.config,self.rows,self.receipt(batch))
        self.assertEqual(result['status'],'nas_verified')
        self.assertTrue((self.sample('S1')/'compact').is_dir())
        self.assertEqual(archive_receipt.matches(self.config,self.rows[0],self.sample('S1'))['destination'],batch['samples'][0]['destination'])

    def test_archive_keeps_qc_metadata_logs_model_bytes_and_is_idempotent(self):
        batch,record=self.archive();root=self.sample('S1')
        self.assertEqual(record['status'],'archived')
        self.assertFalse((root/'compact').exists());self.assertFalse((root/'rna_personal').exists())
        self.assertTrue((root/'qc/final_package_validation.json').exists())
        self.assertTrue((root/'metadata/sample.json').exists());self.assertTrue((root/'logs').is_dir())
        self.assertEqual(json.loads((self.state/'samples/S1.json').read_text())['model_bytes'],8)
        self.assertEqual(transfer.commit(self.config,self.rows,self.receipt(batch),True)['status'],'archived')

    def test_interrupted_source_cleanup_resumes_without_reprocessing(self):
        batch=self.prepare(batch_size=1);remove=shutil.rmtree
        def fail_rna(path):
            if Path(path).name=='rna_personal':raise OSError('의도한 삭제 중단')
            return remove(path)
        with patch.object(transfer.shutil,'rmtree',side_effect=fail_rna),self.assertRaises(OSError):
            transfer.commit(self.config,self.rows,self.receipt(batch),True)
        self.assertTrue(cli.processing_complete(self.sample('S1'),'S1','rna',require_cleanup=True))
        self.assertEqual(archive_receipt.read(self.sample('S1'))['status'],'nas_verified')
        self.assertEqual(transfer.commit(self.config,self.rows,self.receipt(batch),True)['status'],'archived')

    def test_open_source_payload_is_retained_until_reader_closes(self):
        batch=self.prepare(batch_size=1);path=self.sample('S1')/'compact/fake/chr1/0'
        with path.open('rb'),self.assertRaisesRegex(ValueError,'읽는 프로세스'):
            transfer.commit(self.config,self.rows,self.receipt(batch),True)
        self.assertTrue(path.exists())
        self.assertEqual(transfer.commit(self.config,self.rows,self.receipt(batch),True)['status'],'archived')

    def test_archive_is_complete_for_queue_and_worker_after_restart(self):
        self.archive();adapter=local_queue.Adapter(self.config,self.rows,cli)
        self.assertIn('S1',adapter.completed_ids(self.rows))
        with patch('handoff.worker_resources.apply',side_effect=lambda c,e:c),patch.object(cli.subprocess,'run') as run:
            self.assertEqual(cli.worker(self.config,self.rows[0]),0)
            run.assert_not_called()

    def test_archived_sample_is_not_an_explicit_recovery_candidate(self):
        self.archive()
        from handoff.queue_recovery import enqueue
        with self.assertRaisesRegex(ValueError,'NAS 이관 완료'):
            enqueue(self.config,self.rows,['S1'],'01',False,False,transfer.save)

    def test_other_run_or_changed_retained_qc_is_not_silently_accepted(self):
        self.archive()
        with self.assertRaises(ValueError):
            archive_receipt.matches({**self.config,'run_id':'other'},self.rows[0],self.sample('S1'))
        (self.sample('S1')/'qc/final_package_validation.json').write_text('{}')
        with self.assertRaises(ValueError):archive_receipt.read(self.sample('S1'))

    def test_catalog_records_nas_paths_without_claiming_local_availability(self):
        _,record=self.archive();row=archive_receipt.catalog_row(self.sample('S1'))
        self.assertEqual(row['processing_status'],'archived');self.assertEqual(row['data_location'],'nas')
        self.assertEqual(row['qc_all_paths_exist'],'False')
        self.assertEqual(row['compact_root'],record['destination']+'/compact/fake')

    def test_receiver_detects_corruption_and_unexpected_files(self):
        batch=self.prepare(batch_size=1);item=batch['samples'][0];target=self.base/'received'
        shutil.copytree(item['source_root'],target)
        receiver.verify(target,Path(item['inventory_path']),item)
        path=target/'compact/fake/chr1/0';stat=path.stat();path.write_bytes(b'TGCA');os.utime(path,ns=(stat.st_atime_ns,stat.st_mtime_ns))
        with self.assertRaises(ValueError):receiver.verify(target,Path(item['inventory_path']),item)
        path.write_bytes(b'ACGT');(target/'extra').write_text('bad')
        with self.assertRaises(ValueError):receiver.verify(target,Path(item['inventory_path']),item)

    def test_source_and_receiver_reject_symlink_escape(self):
        (self.sample('S1')/'compact/link').symlink_to(self.base/'outside')
        with self.assertRaises(ValueError):self.prepare(batch_size=1)

    def test_watch_requests_ready_samples_without_waiting_for_ten(self):
        args=SimpleNamespace(watch=True,allow_partial=False,batch_size=10,batch=None,
            destination=str(self.nas),source='test-source',ssh=['ssh'],source_code=str(ROOT),
            source_config=str(self.state/'effective_config.json'),source_python=sys.executable)
        response=SimpleNamespace(stdout=json.dumps({'status':'waiting','ready':0,'required':10}))
        with patch.object(receiver.subprocess,'run',return_value=response) as command,patch.object(receiver,'receive') as receive:
            self.assertFalse(receiver.pull_once(args))
            self.assertIn('--allow-partial',command.call_args.args[0][-1])
            receive.assert_not_called()

    def test_watch_waits_only_when_idle_and_immediately_drains_new_completions(self):
        args=SimpleNamespace(watch=True,poll_interval=300)
        with patch.object(receiver,'pull_once',side_effect=[False,True,True,KeyboardInterrupt]) as poll,patch.object(receiver.time,'sleep') as sleep:
            with self.assertRaises(KeyboardInterrupt):receiver.run(args)
            self.assertEqual(poll.call_count,4);sleep.assert_called_once_with(300)

    def test_watch_retries_connection_failure_after_poll_interval(self):
        args=SimpleNamespace(watch=True,poll_interval=300)
        failure=subprocess.CalledProcessError(255,['ssh'],stderr='test connection failure')
        with patch.object(receiver,'pull_once',side_effect=[failure,True,KeyboardInterrupt]) as poll,patch.object(receiver.time,'sleep') as sleep:
            with self.assertRaises(KeyboardInterrupt):receiver.run(args)
            self.assertEqual(poll.call_count,3);sleep.assert_called_once_with(300)

    def test_integrity_errors_stop_watch_and_manual_transfer_does_not_retry(self):
        for watch,error in [(True,ValueError('NAS checksum mismatch')),
                            (False,subprocess.CalledProcessError(255,['ssh']))]:
            with self.subTest(watch=watch),patch.object(receiver,'pull_once',side_effect=error),patch.object(receiver.time,'sleep') as sleep:
                with self.assertRaises(type(error)):receiver.run(SimpleNamespace(watch=watch,poll_interval=300))
                sleep.assert_not_called()

    def test_watch_rejects_pinned_batch_or_invalid_polling_before_remote_call(self):
        common=['pull_completed.py','--source','test-source','--identity','unused',
            '--source-code',str(ROOT),'--source-config',str(self.state/'effective_config.json'),
            '--source-python',sys.executable,'--destination',str(self.nas),'--watch']
        for extra in [['--batch','old_batch'],['--poll-interval','0'],['--batch-size','0']]:
            with self.subTest(extra=extra),patch.object(sys,'argv',common+extra),patch.object(sys,'stderr'),patch.object(receiver.subprocess,'run') as remote:
                with self.assertRaises(SystemExit) as caught:receiver.main()
                self.assertEqual(caught.exception.code,2);remote.assert_not_called()
                self.assertFalse(self.nas.exists())

    def test_second_receiver_cannot_poll_same_nas_destination(self):
        import fcntl
        self.nas.mkdir()
        argv=['pull_completed.py','--source','test-source','--identity','unused',
            '--source-code',str(ROOT),'--source-config',str(self.state/'effective_config.json'),
            '--source-python',sys.executable,'--destination',str(self.nas),'--watch']
        with (self.nas/'.nas_pull.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            with patch.object(sys,'argv',argv),patch.object(sys,'stderr'),patch.object(receiver.subprocess,'run') as remote:
                with self.assertRaises(SystemExit) as caught:receiver.main()
                self.assertEqual(caught.exception.code,2);remote.assert_not_called()
        self.assertTrue((self.sample('S1')/'compact').is_dir())

    def test_full_receiver_with_real_rsync_and_fake_ssh(self):
        batch=self.prepare(batch_size=1);self.nas.mkdir()
        ssh=self.base/'fake_ssh.py';ssh.write_text("import os,sys\nos.execv('/bin/sh',['sh','-c',' '.join(sys.argv[2:])])\n")
        args=SimpleNamespace(destination=str(self.nas),source='local-test',ssh=[sys.executable,str(ssh)],
            source_code=str(ROOT),source_config=str(self.state/'effective_config.json'),
            source_python=sys.executable,remove_source=True,watch=True)
        receiver.receive(args,batch)
        final=Path(batch['samples'][0]['destination'])
        self.assertEqual((final/'compact/fake/chr1/0').read_bytes(),b'ACGT')
        self.assertFalse((self.sample('S1')/'compact').exists())
        self.assertTrue(list((self.nas/'_transfer_records').rglob('*.verified.json')))
        with patch.object(receiver.subprocess,'run',wraps=receiver.subprocess.run) as commands:
            receiver.receive(args,self.prepare(batch_id=batch['batch_id']))
        payload=args.source+':'+batch['samples'][0]['source_root']+'/'
        self.assertFalse(any(payload in call.args[0] for call in commands.call_args_list))
        self.assertEqual(self.prepare()['ready'],1)

    def test_receiver_failure_keeps_source_and_resumes_same_batch(self):
        batch=self.prepare(batch_size=1);self.nas.mkdir()
        ssh=self.base/'fake_ssh.py';ssh.write_text("import os,sys\nos.execv('/bin/sh',['sh','-c',' '.join(sys.argv[2:])])\n")
        args=SimpleNamespace(destination=str(self.nas),source='local-test',ssh=[sys.executable,str(ssh)],
            source_code=str(ROOT),source_config=str(self.state/'effective_config.json'),
            source_python=sys.executable,remove_source=True)
        def corrupt(directory,*unused):
            path=Path(directory)/'compact/fake/chr1/0';stat=path.stat();path.write_bytes(b'TGCA')
            os.utime(path,ns=(stat.st_atime_ns,stat.st_mtime_ns))
            raise ValueError('의도한 수신 검증 실패')
        with patch.object(receiver,'verify',side_effect=corrupt),self.assertRaises(ValueError):
            receiver.receive(args,batch)
        self.assertTrue((self.sample('S1')/'compact').is_dir())
        self.assertFalse((self.sample('S1')/archive_receipt.RECEIPT).exists())
        self.assertFalse(Path(batch['samples'][0]['destination']).exists())
        receiver.receive(args,self.prepare())
        self.assertEqual(archive_receipt.read(self.sample('S1'))['status'],'archived')


if __name__=='__main__':unittest.main()
