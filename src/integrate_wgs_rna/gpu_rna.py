"""기존 STAR 호출 인터페이스를 GPU STAR로 연결하고 실효 옵션 차이를 기록함."""
from pathlib import Path
import hashlib,json,math,os,shutil,subprocess,sys
from .gpu_pool import reserved,GIB
from .parabricks_runtime import Runtime

# 공식 Parabricks4.8 옵션과 현재 STAR 인자의 의미 대응임.
OPTIONS={'--twopassMode':'--two-pass-mode','--outSAMattributes':'--out-sam-attributes',
 '--outSAMunmapped':'--out-sam-unmapped','--outSAMstrandField':'--out-sam-strand-field',
 '--outFilterMultimapNmax':'--max-out-filter-multimap','--outFilterMismatchNmax':'--max-out-filter-mismatch',
 '--outFilterMismatchNoverLmax':'--max-out-filter-mismatch-ratio',
 '--outFilterScoreMinOverLread':'--min-score-filter','--outFilterMatchNminOverLread':'--min-match-filter-normalized',
 '--outFilterIntronMotifs':'--out-filter-intron-motifs',
 '--alignIntronMin':'--min-intron-size','--alignIntronMax':'--max-intron-size',
 '--alignMatesGapMax':'--max-align-mates-gap','--alignSJoverhangMin':'--min-align-sj-overhang',
 '--alignSJDBoverhangMin':'--min-align-sjdb-overhang','--quantMode':'--quantMode'}
FIXED_DEFAULTS={'--outSAMtype':['BAM','SortedByCoordinate'],
                '--alignSoftClipAtReferenceEnds':['Yes'],'--genomeLoad':['NoSharedMemory'],
                '--scoreGenomicLengthLog2scale':['-0.25']}
DIFFERENCES={'--outFilterType':{'cpu':['BySJout'],'gpu':['Normal'],
    'reason':'Parabricks4.8에 BySJout 설정 인터페이스가 없어 GPU STAR 기본값을 기록함'},
    '--limitSjdbInsertNsj':{'cpu':['1200000'],'gpu':['1000000'],
    'reason':'pbrun에는 삽입 SJ 수 제한 인자가 없음. max-collapsed-junctions는 limitOutSJcollapsed이므로 대체하지 않음'}}

def parse(args):
    result={};key=None
    for value in args:
        if value.startswith('--'):
            key=value
            if key in result:raise ValueError('중복 STAR 인자: '+key)
            result[key]=[]
        elif key is None:raise ValueError('STAR 인자 해석 실패: '+value)
        else:result[key].append(value)
    return result

def validate_policy():
    policy=os.environ.get('MFM_GPU_RNA_POLICY','strict')
    if policy not in {'strict','experimental'}:raise ValueError('GPU RNA policy는 strict|experimental이어야 함')
    if policy=='strict':raise RuntimeError('현재 GTEx BySJout은 GPU STAR에서 설정 불가함. CPU를 사용하거나 --rna-gpu-policy experimental로 차이를 명시해야 함')


def effective_contract():
    return {'backend':'parabricks','compatible_star':'2.7.2a','cpu_profile_star':'2.7.11b',
            'profile_equivalent':False,'differences':DIFFERENCES,
            'compatibility_policy':os.environ.get('MFM_GPU_RNA_POLICY','strict'),
            'independent_source_trace':'기존 검증 유지함'}

def index_star():
    path=os.environ.get('MFM_GPU_RNA_STAR',str(Path(__file__).resolve().parents[2]/'tools/gpu/STAR_2.7.2a'))
    if not Path(path).is_file():raise RuntimeError('GPU RNA index용 STAR2.7.2a 없음. 00_environment.sh install-gpu-tools 또는 --gpu-rna-star 필요함')
    version=subprocess.run([path,'--version'],text=True,capture_output=True,check=True).stdout.strip()
    if version!='2.7.2a':raise ValueError('GPU RNA index는 STAR2.7.2a가 필요함: '+version)
    return path

def integrated_markdup_enabled():
    from .parabricks_bam import enabled
    mode=os.environ.get('RNA_GPU_MARKDUP_MODE','separate')
    if mode not in {'integrated','separate'}:raise ValueError('RNA_GPU_MARKDUP_MODE는 integrated|separate이어야 함')
    return os.environ.get('MFM_RNA_BACKEND')=='parabricks' and enabled() and mode=='integrated'


def integrated_markdup_record(directory,output):
    from rna_cache import cache_fields,file_fingerprint
    metrics=Path(directory)/'MarkDuplicates.metrics.txt'
    if not metrics.is_file() or not metrics.stat().st_size:raise ValueError('GPU 통합 duplicate metrics 누락')
    lines=[line for line in metrics.read_text().splitlines() if line.strip() and not line.startswith('#')]
    required={'READ_PAIRS_EXAMINED','READ_PAIR_DUPLICATES','READ_PAIR_OPTICAL_DUPLICATES',
              'UNPAIRED_READS_EXAMINED','UNPAIRED_READ_DUPLICATES','UNMAPPED_READS','PERCENT_DUPLICATION'}
    if len(lines)<2 or not required.issubset(lines[0].split('\t')):raise ValueError('GPU 통합 duplicate metrics 표/필수 열 누락')
    values=dict(zip(lines[0].split('\t'),lines[1].split('\t')))
    if any(float(values[name])<0 for name in required):raise ValueError('GPU 통합 duplicate metrics 값 오류')
    return {'schema':'mfm_rna_integrated_markdup_v1','bam':str(Path(output).resolve()),
            'metrics':str(metrics.resolve()),'remove_duplicates':False,
            'baseline_sort_order':'coordinate','optical_duplicate_pixel_distance':100,
            **cache_fields({'adapter':file_fingerprint(__file__,content=True)},[Path(output),metrics])}


def validate_integrated_markdup(record,bam):
    from rna_cache import cache_valid,file_fingerprint
    if (record.get('schema')!='mfm_rna_integrated_markdup_v1' or
        record.get('remove_duplicates') is not False or
        record.get('baseline_sort_order')!='coordinate' or
        record.get('optical_duplicate_pixel_distance')!=100 or
        Path(record.get('bam','')).resolve()!=Path(bam).resolve() or
        record.get('cache_inputs',{}).get('adapter')!=file_fingerprint(__file__,content=True) or
        not cache_valid(record) or
        [row.get('path') for row in record.get('cache_outputs',[])]!=
            [str(Path(bam).resolve()),str(Path(record.get('metrics','')).resolve())]):
        raise ValueError('GPU 통합 중복표시 계약/입력 변경: 정렬 단계부터 다시 실행 필요함')
    return Path(record['metrics'])


def align_arguments(values,reference,output,threads,sjdb_overhang=None,*,integrated_markdup=False):
    validate_policy()
    required={'--genomeDir','--readFilesIn','--outFileNamePrefix','--runMode','--runThreadN',
              '--readFilesCommand','--outSAMattrRGline'}
    for key,expected in FIXED_DEFAULTS.items():
        if key in values and values[key]!=expected:raise ValueError('GPU STAR 고정 인자 불일치: '+key)
    unsupported=set(values)-required-set(OPTIONS)-set(FIXED_DEFAULTS)-set(DIFFERENCES)
    if unsupported:raise ValueError('GPU STAR 미지원 인자를 변경하지 않음: '+','.join(sorted(unsupported)))
    for key,contract in DIFFERENCES.items():
        if values.get(key,contract['cpu'])!=contract['cpu']:raise ValueError('검토되지 않은 RNA filtering 옵션: '+key)
    fastqs=values['--readFilesIn']
    if len(fastqs)!=2:raise ValueError('GPU RNA는 paired FASTQ2개를 요구함')
    command=['rna_fq2bam','--ref',str(reference),'--genome-lib-dir',values['--genomeDir'][0],
             '--in-fq',*fastqs,'--output-dir',str(output.parent),'--out-bam',str(output),
             '--read-files-command',values.get('--readFilesCommand',['cat'])[0],
             '--num-threads',str(threads),'--num-streams-per-gpu','1','--num-gpus','1',
             '--gpusort','--gpuwrite',
             '--memory-limit',str(int(os.environ.get('MFM_GPU_HOST_MEMORY_BYTES',str(128*GIB)))//10**9)]
    if integrated_markdup:
        command += ['--out-duplicate-metrics',str(output.parent/'MarkDuplicates.metrics.txt'),
                    '--optical-duplicate-pixel-distance','100']
    else:command += ['--no-markdups']
    if sjdb_overhang is not None:command += ['--sjdb-overhang',str(int(sjdb_overhang))]
    for key,flag in OPTIONS.items():
        if key in values:command += [flag,*values[key]]
    rg=dict(v.split(':',1) for v in values.get('--outSAMattrRGline',[]) if ':' in v)
    if 'SM' in rg:command += ['--read-group-sm',rg['SM']]
    if 'ID' in rg:command += ['--read-group-id-prefix',rg['ID']]
    return command

def main(args=None):
    args=list(sys.argv[1:] if args is None else args)
    if args==['--version']:
        print('2.7.2a-parabricks4.8-'+hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:12]);return 0
    values=parse(args);validate_policy();mode=values.get('--runMode',['alignReads'])[0]
    if mode=='genomeGenerate':
        star=index_star();directory=Path(values['--genomeDir'][0]);directory.mkdir(parents=True,exist_ok=True)
        ref=directory/'gpu_reference.fa'
        # 흡수용 contig까지 같은 FASTA에 포함하여 GPU 도구의 --ref와 STAR index를 일치시킴.
        with ref.open('wb') as out:
            for path in values['--genomeFastaFiles']:
                with open(path,'rb') as source:shutil.copyfileobj(source,out,16<<20)
        import pysam
        pysam.faidx(str(ref))
        replaced=[]
        for key,items in values.items():replaced += [key,*( [str(ref)] if key=='--genomeFastaFiles' else items)]
        subprocess.run([star,*replaced],check=True)
        (directory/'gpu_index_contract.json').write_text(json.dumps({**effective_contract(),
            'reference':str(ref.resolve()),'index_star':star,
            'sjdb_overhang':int(values.get('--sjdbOverhang',['100'])[0])},ensure_ascii=False,indent=2)+'\n')
        return 0
    if mode!='alignReads':raise ValueError('지원하지 않는 GPU STAR 실행모드: '+mode)
    runtime=Runtime();runtime.validate()
    index=Path(values['--genomeDir'][0]);contract_path=index/'gpu_index_contract.json'
    if not contract_path.is_file():raise ValueError('GPU STAR2.7.2a index 계약이 없음. RNA index 단계부터 실행 필요함')
    contract=json.loads(contract_path.read_text());reference=Path(contract['reference'])
    directory=Path(values['--outFileNamePrefix'][0]);directory.mkdir(parents=True,exist_ok=True)
    output=directory/'Aligned.sortedByCoord.out.bam';threads=int(values['--runThreadN'][0])
    os.environ['MFM_GPU_HOST_CPU_LIMIT']=str(threads)
    integrated=integrated_markdup_enabled()
    command=align_arguments(values,reference,output,threads,contract['sjdb_overhang'],integrated_markdup=integrated)
    from .gpu_stage_resources import StageBudget
    os.environ.setdefault('MFM_GPU_HAP_ID','1')
    os.environ.setdefault('MFM_GPU_STAGE_CPUS',str(threads))
    budget=StageBudget()
    with budget.gpu(int(os.environ.get('MFM_GPU_ALIGNMENT_VRAM_BYTES',str(24*GIB))),purpose='rna_alignment') as gpu:
        info=runtime.run(command,gpu,[reference,index,directory,*values['--readFilesIn']],log=directory/'gpu_runtime.log')
    budget.wait('gpu_complete')
    for name in ['Aligned.sortedByCoord.out.bam','Log.final.out','ReadsPerGene.out.tab']:
        path=directory/name
        if not path.is_file() or not path.stat().st_size:raise RuntimeError('GPU RNA 필수 결과 누락: '+name)
    markdup=None
    if integrated:
        import pysam
        from .parabricks_bam import normalize_markdup_header
        pysam.quickcheck(str(output))
        normalize_markdup_header(output,directory,directory/'header.log')
        markdup=integrated_markdup_record(directory,output)
    (directory/'gpu_backend.json').write_text(json.dumps({**effective_contract(),**info,'host_resource_phases':budget.events,
        'integrated_markdup':markdup,
        'requested_star_args':values,'effective_index_sjdb_overhang':contract['sjdb_overhang'],
        'effective_parabricks_args':command},ensure_ascii=False,indent=2)+'\n')
    return 0

if __name__=='__main__':raise SystemExit(main())
