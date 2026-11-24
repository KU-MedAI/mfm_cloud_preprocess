"""최종 검증·ASE·출력 경로와 catalog 상태를 일치시킨다."""
import csv, importlib.util, json, subprocess, sys, tempfile, unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('catalog',ROOT/('27_make_catalog.py' if (ROOT/'09_personal_fasta.sh').is_file() else '25_make_catalog.py'))
catalog=importlib.util.module_from_spec(spec);spec.loader.exec_module(catalog)

def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(data) if isinstance(data,dict) else data)

def package(base,ready=True):
    root=base/'S';rna=root/'rna_personal/panel'
    for h in ['HAP1','HAP2']:(rna/h/'coverage.zarr').mkdir(parents=True)
    for name in ['allelic_counts.tsv','ase_training_eligibility.tsv']:write(rna/name,'x\n')
    write(rna/'ase_qc.json',{'sample_qc':{'sample_training_ready':ready}})
    write(root/'metadata/sample.json',{'sequence_variant_policy':'snv_indel','dataset_role':'post_training','layout_name':'layout','rna_personal':{'coverage_zarr':{h:str(rna/h/'coverage.zarr') for h in ['HAP1','HAP2']},'allelic_counts':str(rna/'allelic_counts.tsv'),'ase_qc_json':str(rna/'ase_qc.json'),'ase_training_eligibility':str(rna/'ase_training_eligibility.tsv')}})
    write(root/'compact/layout/manifest.json',{'status':'complete','validation_level':'full','records':{c:{'contig':c,'status':'pass'} for c in catalog.CHROMS}})
    write(root/'qc/stage12_trace/all_to_both/trace_manifest.json',{'status':'pass','compact_root':str(root/'compact/layout'),'contigs':{c:'pass' for c in catalog.CHROMS}})
    write(root/'qc/final_package_validation.json',{'status':'pass','failed_checks':[]})
    return root,rna

class CatalogStatus(unittest.TestCase):
    def test_invalid_final_reports_are_not_eligible(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,_=package(Path(tmp))
            for report in [{'status':'fail','failed_checks':[]},{'status':'pass','failed_checks':['error']},{'status':'pass','failed_checks':['']},{'status':'pass','failed_checks':False},{'status':'pass'},{}]:
                write(root/'qc/final_package_validation.json',report);row=catalog.scan_sample(root,'TEST',{}, {})
                self.assertEqual(row['training_eligible'],'False');self.assertNotEqual(row['processing_status'],'complete')
    def test_ase_true_false_and_types(self):
        for value,expected in [(True,'True'),(False,'False'),('true','False'),(1,'False')]:
            with tempfile.TemporaryDirectory() as tmp:
                root,_=package(Path(tmp),value);row=catalog.scan_sample(root,'TEST',{}, {})
                self.assertEqual(row['training_eligible'],expected)
                self.assertEqual(row['processing_status'],'complete' if expected=='True' else 'processed_qc_fail')
    def test_validation_pending_for_missing_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,_=package(Path(tmp));(root/'qc/final_package_validation.json').unlink()
            self.assertEqual(catalog.scan_sample(root,'TEST',{}, {})['processing_status'],'validation_pending')
    def test_final_output_excludes_missing_model_and_failed_validation(self):
        for case in ['failed','missing_model']:
            with tempfile.TemporaryDirectory() as tmp:
                base=Path(tmp);root,rna=package(base)
                if case=='failed':write(root/'qc/final_package_validation.json',{'status':'fail','failed_checks':['error']})
                else:(rna/'allelic_counts.tsv').unlink()
                out=base/'catalog.tsv';r=subprocess.run([sys.executable,str(ROOT/('27_make_catalog.py' if (ROOT/'09_personal_fasta.sh').is_file() else '25_make_catalog.py')),str(base),'--out',str(out)],capture_output=True,text=True)
                self.assertEqual(r.returncode,1,r.stderr)
                with out.open() as f:self.assertEqual(list(csv.DictReader(f,delimiter='\t')),[])
    def test_ase_qc_fail_is_retained_as_ineligible_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp);package(base,False);out=base/'catalog.tsv'
            r=subprocess.run([sys.executable,str(ROOT/('27_make_catalog.py' if (ROOT/'09_personal_fasta.sh').is_file() else '25_make_catalog.py')),str(base),'--out',str(out)],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
            with out.open() as f:row=next(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(row['training_eligible'],'False');self.assertEqual(row['processing_status'],'processed_qc_fail')
            self.assertIn('[발견한 전체 샘플]',r.stdout);self.assertIn('[최종 catalog 출력 샘플]',r.stdout)

if __name__=='__main__':unittest.main()
