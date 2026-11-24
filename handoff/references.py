"""참조 역할을 준비하고 CRAM @SQ 길이/M5에 맞는 FASTA를 선택한다."""
from __future__ import annotations
import gzip,hashlib,json,os,subprocess
from pathlib import Path
from .downloads import fetch

NCBI='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
IGSR='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa'
BROAD='https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta'
GENCODE='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz'
HGSVC='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz'


def needs_hgsvc_reference(row):
    """HGSVC PAV callset은 NCBI analysis set과 N/IUPAC 마스킹이 다른 참조를 쓴다."""
    return row.get('vcf_source') == 'public' and row.get('project_id') in {'HGSVC2', 'HGSVC3'}


def variant_reference(row, refs):
    if needs_hgsvc_reference(row):
        if not refs.get('hgsvc'):
            raise ValueError('HGSVC public VCF에는 전용 hgsvc 참조가 필요합니다.')
        return refs['hgsvc']
    return refs['autosomes']

def index(path,samtools):
    p=Path(path)
    if not Path(str(p)+'.fai').exists():subprocess.run([samtools,'faidx',str(p)],check=True)
    return p

def decompress(source,target):
    target=Path(target)
    if target.exists():return target
    tmp=Path(str(target)+'.part')
    with gzip.open(source,'rb') as src,tmp.open('wb') as dst:
        import shutil;shutil.copyfileobj(src,dst,8<<20)
    os.replace(tmp,target);mark_derived(target);return target

def mark_derived(target):
    target=Path(target)
    Path(str(target)+'.owned.json').write_text(json.dumps({'schema':'mfm_derived_reference_v1','path':str(target.resolve())})+'\n')

def write_subset(source,target,names):
    import pysam
    target=Path(target);tmp=Path(str(target)+'.part')
    with pysam.FastaFile(str(source)) as src,tmp.open('w') as dst:
        for name in names:
            sequence=src.fetch(name);dst.write('>'+name+'\n')
            for i in range(0,len(sequence),60):dst.write(sequence[i:i+60]+'\n')
    os.replace(tmp,target);mark_derived(target);return target

def prepare(config,root,needs_rna=False,needs_illumina=False,needs_hgsvc=False):
    import fcntl,pysam
    from .bundled_assets import references as bundled_references
    bundled = bundled_references(config, root, needs_rna, needs_illumina, needs_hgsvc)
    if bundled is not None:
        return bundled
    root=Path(root);data=Path(config['data_root']);state=Path(config['state_root'])
    ref=data/'references';ref.mkdir(parents=True,exist_ok=True)
    samtools=config.get('samtools','samtools')
    hgsvc=None
    with (ref/'prepare.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        ncbi=ref/'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
        if not ncbi.exists():
            gz=fetch({'url':config.get('ncbi_reference_url',NCBI),'path':'references/'+ncbi.name+'.gz',
                'size_bytes':872949833,'checksum_type':'md5','checksum':'a08035b6a6e31780e96a34008ff21bd6'},data,state)
            decompress(gz,ncbi)
        index(ncbi,samtools)
        autosomes=ref/'hg38_autosomes.fa'
        if not autosomes.exists():write_subset(ncbi,autosomes,['chr'+str(i) for i in range(1,23)])
        index(autosomes,samtools)
        if needs_hgsvc:
            hgsvc=ref/'hgsvc_hg38_autosomes.fa'
            if not hgsvc.exists():
                # 원본은 공개 체크섬으로 검증하고, 역할이 다른 기존 참조를 덮어쓰지 않는다.
                source=ref/'hgsvc_hg38.no_alt.fa'
                if not source.exists():
                    gz=fetch({'url':HGSVC,'path':'references/'+source.name+'.gz',
                        'size_bytes':963047116,'checksum_type':'md5','checksum':'6a5ceceed91ae5c14026454f09c4b6c4',
                        'shared':True,'kind':'reference'},data,state)
                    decompress(gz,source)
                index(source,samtools)
                write_subset(source,hgsvc,['chr'+str(i) for i in range(1,23)])
            index(hgsvc,samtools)
        broad=ref/'Homo_sapiens_assembly38.fasta'
        if needs_rna or needs_illumina:
            if not broad.exists():fetch({'url':config.get('broad_reference_url',BROAD),'path':'references/'+broad.name,
                'size_bytes':3249912778,'checksum_type':'md5','checksum':'7ff134953dcca8c8997453bbb80b6b5e'},data,state)
            index(broad,samtools)
        if needs_illumina:
            igsr=ref/'GRCh38_full_analysis_set_plus_decoy_hla.fa'
            if not igsr.exists():fetch({'url':config.get('illumina_reference_url',IGSR),'path':'references/'+igsr.name},data,state)
            index(igsr,samtools)
        if needs_rna:
            gtf=ref/'gencode.v47.annotation.gtf.gz'
            if not gtf.exists():fetch({'url':config.get('gencode_url',GENCODE),'path':'references/'+gtf.name},data,state)
            extra=ref/'rna_index_nonautosome_gtex.fa'
            if not extra.exists():
                with pysam.FastaFile(str(broad)) as fasta:
                    names=[n for n in fasta.references if n not in {'chr'+str(i) for i in range(1,23)}
                      and not n.endswith(('_alt','_decoy')) and not n.startswith('HLA')]
                write_subset(broad,extra,names)
            index(extra,samtools)
    result={'ncbi':str(ncbi),'autosomes':str(autosomes),'broad':str(broad),
        'illumina':str(ref/'GRCh38_full_analysis_set_plus_decoy_hla.fa'),
        'gtf':str(ref/'gencode.v47.annotation.gtf.gz'),'rna_extra':str(ref/'rna_index_nonautosome_gtex.fa')}
    if hgsvc is not None:result['hgsvc']=str(hgsvc)
    return result

def exact_cram_reference(cram,candidates,out,samtools='samtools'):
    import fcntl,pysam
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    header=subprocess.check_output([samtools,'view','-H',str(cram)],text=True)
    sq=[]
    for line in header.splitlines():
        if line.startswith('@SQ\t'):
            tags={s[:2]:s[3:] for s in line.split('\t')[1:] if s[2:3]==':'}
            if not tags.get('M5'):raise ValueError('CRAM @SQ M5 누락: '+tags.get('SN',''))
            sq.append(tags)
    fingerprint=hashlib.sha256(json.dumps(sq,sort_keys=True).encode()).hexdigest()
    out=out.with_name(out.stem+'_'+fingerprint[:20]+'.fa');sidecar=Path(str(out)+'.json')
    with Path(str(out)+'.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if out.exists() and sidecar.exists() and json.loads(sidecar.read_text()).get('sq_fingerprint')==fingerprint:
            return index(out,samtools)
        handles=[pysam.FastaFile(str(p)) for p in candidates if Path(p).is_file()]
        tmp=Path(str(out)+'.part');selected=[]
        try:
            with tmp.open('w') as dst:
                for tags in sq:
                    name=tags['SN'];found=False
                    for ref in handles:
                        if name not in ref.references or ref.get_reference_length(name)!=int(tags['LN']):continue
                        seq=ref.fetch(name).upper()
                        if hashlib.md5(seq.encode()).hexdigest()!=tags['M5'].lower():continue
                        dst.write('>'+name+'\n')
                        for i in range(0,len(seq),60):dst.write(seq[i:i+60]+'\n')
                        selected.append({'contig':name,'length':len(seq),'m5':tags['M5'],'source':ref.filename.decode() if isinstance(ref.filename,bytes) else ref.filename})
                        found=True;break
                    if not found:raise ValueError('CRAM M5 일치 reference 없음: '+name+' M5='+tags['M5']+'; config의 extra_decode_references에 검증된 FASTA를 추가하세요.')
            os.replace(tmp,out);mark_derived(out)
            sidecar.write_text(json.dumps({'sq_fingerprint':fingerprint,'contigs':selected,'status':'pass'},ensure_ascii=False,indent=2)+'\n')
        finally:
            for ref in handles:ref.close()
            if tmp.exists():tmp.unlink()
    return index(out,samtools)
