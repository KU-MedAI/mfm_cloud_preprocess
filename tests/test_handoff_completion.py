"""original 중간 산출물 보존과 전달본 cleanup 완료 판정의 차이를 검증한다."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from handoff import cli

class CompletionPolicy(unittest.TestCase):
    def test_original_requires_final_pass_and_preserved_assets(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(cli,'ORIGINAL_LAYOUT',True):
            root=Path(tmp)/'S';(root/'qc').mkdir(parents=True)
            report={'schema':'mfm_final_package_validation_v1','sample':'S','mode':'rna',
                    'status':'pass','failed_checks':[]}
            target=root/'qc/final_package_validation.json';target.write_text(json.dumps(report))
            self.assertFalse(cli.processing_complete(root,'S','rna'))
            (root/'individual_fasta').mkdir();(root/'individual_fasta/permanent_assets.json').write_text('{}')
            (root/'realigned').mkdir()
            self.assertTrue(cli.processing_complete(root,'S','rna'))
            self.assertFalse(cli.processing_complete(root,'S','rna',require_cleanup=True))
            for change in ({'sample':'OTHER'},{'mode':'wgs'},{'status':'fail'},
                           {'failed_checks':['hash_mismatch']},{'failed_checks':None}):
                target.write_text(json.dumps({**report,**change}))
                self.assertFalse(cli.processing_complete(root,'S','rna'))
            target.write_text('{broken')
            self.assertFalse(cli.processing_complete(root,'S','rna'))

    def test_delivery_success_still_requires_completed_cleanup(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(cli,'ORIGINAL_LAYOUT',False):
            root=Path(tmp)/'S';(root/'qc').mkdir(parents=True)
            (root/'qc/final_package_validation.json').write_text(json.dumps({'status':'pass'}))
            self.assertFalse(cli.processing_complete(root,'S','rna'))
            meta=root/'metadata';meta.mkdir()
            completed={'schema':'mfm_delivery_cleanup_v1','sample':'S','mode':'rna','dry_run':False,
                       'status':'complete','post_cleanup_validation':{
                           'schema':'mfm_delivery_final_validation_v1','sample':'S','mode':'rna',
                           'status':'pass','errors':[]}}
            (meta/'cleanup_completed.json').write_text(json.dumps(completed))
            self.assertTrue(cli.processing_complete(root,'S','rna'))
            (meta/'cleanup_in_progress.json').write_text('{}')
            self.assertFalse(cli.processing_complete(root,'S','rna'))

    def test_original_model_only_cleanup_report_controls_completion(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(cli,'ORIGINAL_LAYOUT',True):
            root=Path(tmp)/'S';meta=root/'metadata';meta.mkdir(parents=True)
            completed={'schema':'mfm_delivery_cleanup_v1','sample':'S','mode':'rna','dry_run':False,
                       'status':'complete','post_cleanup_validation':{
                           'schema':'mfm_delivery_final_validation_v1','sample':'S','mode':'rna',
                           'status':'pass','errors':[]}}
            (meta/'cleanup_completed.json').write_text(json.dumps(completed))
            self.assertTrue(cli.processing_complete(root,'S','rna'))
            (meta/'cleanup_in_progress.json').write_text('{}')
            self.assertFalse(cli.processing_complete(root,'S','rna'))

    def test_original_model_only_config_allows_owned_raw_cleanup(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(cli,'ORIGINAL_LAYOUT',True):
            cfg=Path(tmp)/'config.json'
            cfg.write_text(json.dumps({'data_root':'input','output_root':'output','state_root':'state',
                'resource_profile':'balanced','local_jobs':1,'sample_cpus':4,'elastic_max_jobs':1,
                'remove_owned_raw_inputs':True,'cleanup_shared_fasta':True,
                'pipeline_env':{'MFM_MODEL_ONLY_CLEANUP':'1','CLEAN_INTERMEDIATES_ON_SUCCESS':'1'},
                'slurm':{'cpus':4,'concurrency':1}}))
            c=cli.load_config(cfg)
            self.assertTrue(c['remove_owned_raw_inputs'])
            self.assertFalse(c['cleanup_shared_fasta'])

    def test_original_config_disables_raw_and_reference_deletion(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(cli,'ORIGINAL_LAYOUT',True):
            root=Path(tmp);cfg=root/'config.json'
            cfg.write_text(json.dumps({'data_root':'input','output_root':'output','state_root':'state',
                'resource_profile':'balanced','local_jobs':1,'sample_cpus':4,'elastic_max_jobs':1,
                'remove_owned_raw_inputs':True,'cleanup_shared_fasta':True,
                'slurm':{'cpus':4,'concurrency':1}}))
            c=cli.load_config(cfg)
            self.assertFalse(c['remove_owned_raw_inputs'])
            self.assertFalse(c['cleanup_shared_fasta'])

if __name__=='__main__':unittest.main()
