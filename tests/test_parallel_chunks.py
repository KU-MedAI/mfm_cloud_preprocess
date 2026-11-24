"""bounded QC·공유CPU차용·실제FASTA순차/병렬동등성을검증함."""
from pathlib import Path
from unittest.mock import patch
import hashlib
import json
import os
import shutil
import tempfile
import unittest

import numpy as np
from numcodecs import blosc
from integrate_wgs_rna.compact.chunk_validation import validate_chunk,validate_cigars,borrowed_threads
from integrate_wgs_rna.direct_compact_writer import _ReadValidationChunkCache as ReadValidationChunkCache
from integrate_wgs_rna.memory_policy import GIB


class ChunkTests(unittest.TestCase):
    def data(self,n=100):
        return {'cigar_offset':np.arange(n,dtype=np.int64)*3,'cigar_len':np.full(n,3,dtype=np.int32),
                'ref_start':np.arange(n,dtype=np.int64)*100,'ref_end':np.arange(n,dtype=np.int64)*100+12,
                'lengths':np.full(n,15,dtype=np.int32),'segment_type':np.zeros(n,dtype=np.int8)},np.tile(np.asarray([(5<<4)|0,(3<<4)|1,(7<<4)|0],dtype=np.uint32),n)

    def test_벡터와기존개별decode동등(self):
        arrays,ops=self.data();cache=ReadValidationChunkCache(arrays,ops,100)
        for i in range(100):cache.decode(i)
        self.assertEqual(validate_cigars(arrays,ops,100)['records'],100)

    def test_오류op_query_reference길이를양쪽모두거부(self):
        for kind in ['op','length','query','reference']:
            a,ops=self.data(3)
            if kind=='op':ops[3]=(5<<4)|4
            if kind=='length':ops[3]=0
            if kind=='query':a['lengths'][1]=999
            if kind=='reference':a['ref_end'][1]=999
            with self.assertRaises(ValueError):ReadValidationChunkCache(a,ops,3).decode(1)
            with self.assertRaises(ValueError):validate_cigars(a,ops,3)

    def test_비연속offset과미참조op도기존계약유지(self):
        a,ops=self.data(3);ops=np.insert(ops,3,np.uint32(15));a['cigar_offset'][1:]+=1
        for i in range(3):ReadValidationChunkCache(a,ops,3).decode(i)
        self.assertEqual(validate_cigars(a,ops,3)['records'],3)

    def test_과대chunk를분할하고레코드누락없음(self):
        a,ops=self.data(100)
        with patch('integrate_wgs_rna.compact.chunk_validation.MAX_OPS',30):r=validate_cigars(a,ops,100)
        self.assertGreater(r['chunks'],1);self.assertEqual(r['records'],100)

    def test_단일거대CIGAR는기존decode로안전검증(self):
        a,ops=self.data(1)
        with patch('integrate_wgs_rna.compact.chunk_validation.MAX_OPS',1):self.assertEqual(validate_cigars(a,ops,1)['records'],1)

    def test_독립호출도작은SlurmRAM과CPU상한준수(self):
        with patch.dict(os.environ,{'MFM_SAMPLE_CPU_BUDGET':'68','SLURM_CPUS_PER_TASK':'2','SLURM_MEM_PER_NODE':'4096'},clear=True):
            with borrowed_threads(100) as workers:self.assertEqual(workers,1)

    def test_추가thread는남는CPU_RAM만예약하고반환(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'leases').mkdir();(root/'waiters').mkdir();owner=root/'leases/owner.json';owner.write_text(json.dumps({'sample':'S','cpus':4,'total_cpus':8,'memory_bytes':8*GIB,'manager_pid':os.getpid(),'child_pid':None}))
            prior=blosc.get_nthreads();blosc.set_nthreads(1)
            try:
                with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(owner),'MFM_RAM_LIMIT_BYTES':str(100*GIB),'MFM_QC_CHUNK_THREADS':'auto'}):
                    with borrowed_threads(20) as workers:
                        self.assertLessEqual(workers,5)
                        rows=[json.loads(p.read_text()) for p in (root/'leases').glob('*.json')]
                        self.assertLessEqual(sum(r['cpus'] for r in rows),8)
                        self.assertTrue(all(r.get('child_pid') is None for r in rows))
                    self.assertEqual(len(list((root/'leases').glob('*.json'))),1)
            finally:blosc.set_nthreads(prior)


class FastaTests(unittest.TestCase):
    def test_실제bcftools_병렬FASTA_chain_QC동일(self):
        import pysam
        from integrate_wgs_rna.personal_fasta import PersonalFastaConfig,run_personal_fasta_generation
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);reference=root/'ref.fa';reference.write_text('>chr1\n'+'ACGT'*100+'\n>chr2\n'+'TGCA'*100+'\n');pysam.faidx(str(reference));seq=root/'vcf';seq.mkdir();vcf=seq/'S.seq.vcf';vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=400>\n##contig=<ID=chr2,length=400>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tC\tCT\t.\tPASS\t.\tGT\t0|1\nchr2\t20\t.\tA\tAC\t.\tPASS\t.\tGT\t1|0\n');pysam.tabix_compress(str(vcf),str(vcf)+'.gz',force=True);pysam.tabix_index(str(vcf)+'.gz',preset='vcf',force=True)
            outputs=[]
            for name,jobs in [('serial','1'),('parallel','auto')]:
                output=root/name;cfg=PersonalFastaConfig(cohort='T',seq_vcf_dir=seq,output_root=output,reference_fasta=reference,chroms=['chr1','chr2'],sequence_variant_policy='snv_indel',samples=['S'],jobs=1,bcftools_bin=shutil.which('bcftools'),samtools_bin=shutil.which('samtools'))
                with patch.dict(os.environ,{'MFM_FASTA_CONTIG_JOBS':jobs,'MFM_SAMPLE_CPU_BUDGET':'8','MFM_RAM_LIMIT_BYTES':str(1000*GIB),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}):run_personal_fasta_generation(cfg)
                outputs.append(output/'S/individual_fasta')
            for hap in ['HAP1','HAP2']:
                for chrom in ['chr1','chr2']:
                    for suffix in ['fa','chain']:
                        self.assertEqual((outputs[0]/hap/f'{chrom}.{suffix}').read_bytes(),(outputs[1]/hap/f'{chrom}.{suffix}').read_bytes())
            a=json.loads((outputs[0]/'personal_fasta_qc.json').read_text());b=json.loads((outputs[1]/'personal_fasta_qc.json').read_text());self.assertEqual(a['length_qc'],b['length_qc']);self.assertEqual(a['overlap_skip_count_total'],b['overlap_skip_count_total'])


if __name__=='__main__':unittest.main()
