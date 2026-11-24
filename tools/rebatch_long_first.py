"""MedAI long/pair 우선 배치와파생목록을갱신함. --write 없으면검토만함."""
from pathlib import Path
from collections import Counter
import argparse,csv,io,json,os,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from handoff.batch_priority import rebatch
from handoff.planner import batch_composition

def read(path):
    with path.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
def write(path,rows,fields):
    path.parent.mkdir(parents=True,exist_ok=True);s=io.StringIO(newline='');w=csv.DictWriter(s,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows);tmp=path.with_suffix(path.suffix+'.part');tmp.write_text(s.getvalue());os.replace(tmp,path)
def apply(root,do_write=False):
    source=root/'manifests/분담';split=root/'manifests/리드별';old=read(source/'전체_분담배치.tsv');rows=rebatch(old);fields=list(old[0]);composition=batch_composition(rows)
    report={'samples':len(rows),'counts':dict(Counter(r['assigned_lab'] for r in rows)),'batch_composition':[r for r in composition if r['assigned_lab']=='MedAI'],'first_batch_long_pairs':[r['sample_id'] for r in rows if r['assigned_lab']=='MedAI' and r['read_type']=='long_read' and r['data_group']=='wgs_rna_pair'],'tgil_unchanged':[r for r in rows if r['assigned_lab']=='TGIL']==[r for r in old if r['assigned_lab']=='TGIL'],'written':do_write}
    assert report['samples']==3212 and report['tgil_unchanged']
    if not do_write:return report
    write(source/'전체_분담배치.tsv',rows,fields)
    write(source/'batch_composition.tsv',composition,list(composition[0]))
    for lab in ['MedAI','TGIL']:
        group=[r for r in rows if r['assigned_lab']==lab];slug=lab.lower();write(source/(slug+'_samples.tsv'),group,fields)
        for b in sorted({r['batch_id'] for r in group}):write(source/slug/('batch_'+b+'.tsv'),[r for r in group if r['batch_id']==b],fields)
        for kind in ['long_read','short_read']:
            write(split/(slug+'_'+kind+'.tsv'),[r for r in group if r['read_type']==kind],fields)
            for b in sorted({r['batch_id'] for r in group}):write(split/slug/b/(kind+'.tsv'),[r for r in group if r['read_type']==kind and r['batch_id']==b],fields)
    for kind in ['long_read','short_read']:write(split/('전체_'+kind+'.tsv'),[r for r in rows if r['read_type']==kind],fields)
    summary=[{k:r[k] for k in ['assigned_lab','batch_id','long_read_count','short_read_count','sample_count']} for r in composition];write(split/'배치별_요약.tsv',summary,list(summary[0]))
    lookup={r['sample_id']:r for r in rows};ram=read(split/'RAM재검토.tsv')
    for r in ram:r['batch_id']=lookup[r['sample_id']]['batch_id']
    write(split/'RAM재검토.tsv',ram,list(ram[0]))
    p=source/'분담요약.json';d=json.loads(p.read_text());d.update(batching_strategy='MedAI:long_pair_first_v1;TGIL:deterministic_stratified_best_fit_v2',batch_assignment_seed='MedAI:mfm-long-pair-first-261027;TGIL:mfm-stratified-batches-v2',priority_note='MedAI long89명은001에배정, 그중pair10명먼저. 기관배정과각배치인원보존함.');p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write',action='store_true');a=p.parse_args();print(json.dumps(apply(ROOT,a.write),ensure_ascii=False,indent=2))
