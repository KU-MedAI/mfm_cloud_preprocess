from __future__ import annotations

import argparse
import json
from pathlib import Path

from .cohort_vcf import CohortVcfConfig, run_cohort_vcf_extract
from .assembly_anchor import AssemblyAnchorConfig, run_assembly_anchor
from .denovo_assembly import DenovoAssemblyConfig, run_denovo_assembly
from .direct_compact_writer import DEFAULT_COMPACT_COMPRESSOR, DEFAULT_COMPACT_COMPRESSION_LEVEL, SUPPORTED_COMPACT_COMPRESSORS
from .full_preprocess import FULL_LAYOUT_NAME, FullManifestConfig, FullPreprocessConfig, run_full_manifest, run_full_preprocess
from .haplotag import HaplotagConfig, run_haplotag_generation
from .igsr_ont_manifest import IgsrOntManifestConfig, build_igsr_ont_manifest
from .input_manifest import validate_wgs_input_manifest, write_contract_json
from .legacy_shared_hap_converter import LegacyConvertConfig, run_legacy_convert
from .longread_bamonly_phase import BamOnlyPhaseConfig, run_bamonly_phase
from .longread_clair3 import Clair3LongReadConfig, run_clair3_longread_generation
from .clair3_qc import Clair3QcConfig, run_clair3_qc
from .manifest_runner import ManifestContigRunConfig, ManifestRunConfig, run_manifest, run_manifest_contigs
from .personal_fasta import PersonalFastaConfig, run_personal_fasta_generation
from .personal_realign import PersonalRealignConfig, run_personal_realign_generation
from .public_vcf import PublicVcfConfig, prepare_public_vcf
from .readset_prototype import PrototypeConfig, run_prototype
from .sample_vcf import PerSampleVcfConfig, run_per_sample_vcf_extract
from .window_index import DEFAULT_BLACKLIST_BED, WindowIndexConfig, run_window_index_generation
from .stage70_trace_validator import validate_stage70_trace_root
from .validate_compact import CompactValidationConfig, run_compact_validation


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Integrated WGS/RNA preprocessing prototype")
    subparsers = parser.add_subparsers(dest="command", required=True)

    input_manifest = subparsers.add_parser(
        "validate-wgs-input-manifest",
        help="Validate legacy alignment or multi-file WGS input contract for one sample",
    )
    input_manifest.add_argument("--manifest", required=True, type=Path)
    input_manifest.add_argument("--sample", required=True)
    input_manifest.add_argument("--read-files-manifest", type=Path)
    input_manifest.add_argument("--allow-missing-files", action="store_true")
    input_manifest.add_argument("--output-json", type=Path)

    proto = subparsers.add_parser("prototype-readset", help="Build compact readset encoding prototype for one contig")
    proto.add_argument("--sample-id", required=True)
    proto.add_argument("--alignment", required=True, type=Path)
    proto.add_argument("--contig", required=True)
    proto.add_argument("--output-dir", required=True, type=Path)
    proto.add_argument("--min-mapq", type=int, default=30)
    proto.add_argument("--max-reads", type=int)
    proto.add_argument("--keep-duplicates", action="store_true")
    proto.add_argument("--exclude-supplementary", action="store_true")

    legacy = subparsers.add_parser(
        "convert-legacy-shared-hap",
        help="Convert existing shared-HAP zarr output into the compact integrated 4bit_acgtnp layout",
    )
    legacy.add_argument("--sample-id", required=True)
    legacy.add_argument("--source-sample-id", required=True)
    legacy.add_argument("--legacy-root", required=True, type=Path)
    legacy.add_argument("--output-sample-root", required=True, type=Path)
    legacy.add_argument("--dataset-role", required=True, choices=["pre_training", "post_training"])
    legacy.add_argument("--read-type", required=True, choices=["long_read", "short_read"])
    legacy.add_argument("--project-id", required=True)
    legacy.add_argument("--layout-name", default="legacy_shared_hap_4bit_acgtnp")
    legacy.add_argument("--contigs", help="Comma-separated contigs to convert. Default: all chr*.zarr under legacy root")
    legacy.add_argument("--overwrite", action="store_true")
    legacy.add_argument("--chunk-records", type=int, default=1024)
    legacy.add_argument("--hap-invalid-policy", choices=["fail", "to_n"], default="to_n")
    legacy.add_argument("--rna-bigwig", action="append", type=Path, default=[])
    legacy.add_argument("--source-wgs", action="append", type=Path, default=[])
    legacy.add_argument("--source-vcf", action="append", type=Path, default=[])

    manifest = subparsers.add_parser(
        "convert-legacy-manifest",
        help="Convert multiple legacy shared-HAP samples listed in a TSV manifest",
    )
    manifest.add_argument("--manifest", required=True, type=Path)
    manifest.add_argument("--status-tsv", required=True, type=Path)
    manifest.add_argument("--overwrite", action="store_true")
    manifest.add_argument("--jobs", type=int, default=1)

    contig_manifest = subparsers.add_parser(
        "convert-legacy-manifest-contigs",
        help="Convert legacy shared-HAP manifest with chromosome-level parallelism",
    )
    contig_manifest.add_argument("--manifest", required=True, type=Path)
    contig_manifest.add_argument("--status-tsv", required=True, type=Path)
    contig_manifest.add_argument("--overwrite", action="store_true")
    contig_manifest.add_argument("--jobs", type=int, default=24)


    cohort_vcf = subparsers.add_parser(
        "extract-cohort-common-vcf",
        help="Extract process_code-style cohort common VCFs from a chromosome panel VCF",
    )
    cohort_vcf.add_argument("--cohort", required=True, help="Cohort/project label used in output file names")
    cohort_vcf.add_argument("--sample-list", type=Path, help="Text file with one sample ID per line")
    cohort_vcf.add_argument("--sample", action="append", default=[], help="Sample ID. May be repeated; used when --sample-list is not provided")
    cohort_vcf.add_argument("--output-dir", required=True, type=Path, help="Directory for <cohort>.<contig>.vcf.gz outputs")
    cohort_vcf.add_argument("--chroms", default="chr1-chr22", help="Chromosomes to extract, e.g. chr22 or chr1-chr22")
    cohort_vcf.add_argument(
        "--panel-base",
        default="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV",
        help="Base URL or directory containing chromosome panel VCFs",
    )
    cohort_vcf.add_argument(
        "--panel-file-template",
        default="1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
        help="Panel file template. Supports {chrom}=1, {contig}=chr1, or legacy %%s",
    )
    cohort_vcf.add_argument("--filter-status", default="PASS,.", help="bcftools -f value. process_code default: PASS,.")
    cohort_vcf.add_argument("--min-ac", type=int, default=1, help="bcftools -c value. process_code default: 1")
    cohort_vcf.add_argument("--jobs", type=int, default=6, help="Chromosome-level parallel jobs")
    cohort_vcf.add_argument("--threads-per-chrom", type=int, default=2, help="bcftools --threads per chromosome")
    cohort_vcf.add_argument("--bcftools-bin", default="bcftools")
    cohort_vcf.add_argument("--summary-tsv", type=Path)
    cohort_vcf.add_argument("--metadata-json", type=Path)
    cohort_vcf.add_argument("--idx-cache-dir", type=Path, help="Directory where bcftools may cache remote index files")
    cohort_vcf.add_argument("--overwrite", action="store_true")

    sample_vcf = subparsers.add_parser(
        "extract-per-sample-vcfs",
        help="Create process_code-style per-sample seq/tag VCFs from cohort common VCFs",
    )
    sample_vcf.add_argument("--cohort", required=True, help="Cohort/project label used in file names")
    sample_vcf.add_argument("--common-vcf-dir", required=True, type=Path, help="Directory containing <cohort>.<contig>.vcf.gz common VCFs")
    sample_vcf.add_argument("--output-dir", required=True, type=Path, help="Directory for <sample>.seq/tag.vcf.gz outputs")
    sample_vcf.add_argument("--sample-list", type=Path, help="Text file with one sample ID per line")
    sample_vcf.add_argument("--sample", action="append", default=[], help="Sample ID. May be repeated; used when --sample-list is not provided")
    sample_vcf.add_argument("--chroms", default="chr1-chr22", help="Chromosomes to process, e.g. chr22 or chr1-chr22")
    sample_vcf.add_argument("--sequence-variant-policy", required=True, choices=["snp_only", "snv_indel"], help="Sequence VCF policy")
    sample_vcf.add_argument("--common-vcf-pattern", default="", help="Optional common VCF file template. Supports {cohort}, {contig}, {chrom}")
    sample_vcf.add_argument("--extraction-mode", choices=["chromosome", "concat"], default="chromosome", help="개인 염색체 선추출(기본) 또는 기존 전체 concat")
    sample_vcf.add_argument("--concat-dir", type=Path, help="Directory for <cohort>.autosomes.vcf.gz. Default: --common-vcf-dir")
    sample_vcf.add_argument("--bcftools-bin", default="bcftools")
    sample_vcf.add_argument("--jobs", type=int, default=8, help="Sample-level parallel jobs")
    sample_vcf.add_argument("--concat-threads", type=int, default=0, help="최종 압축 스레드; 0은 샘플 CPU 예산 내 자동 계산")
    sample_vcf.add_argument("--summary-tsv", type=Path)
    sample_vcf.add_argument("--metadata-json", type=Path)
    sample_vcf.add_argument("--overwrite", action="store_true")
    sample_vcf.add_argument(
        "--fast-validate", action="store_true",
        help=("Smoke-only: skip sequence/indel scans; the complete phased tag-GT "
              "contract is still checked for every record"),
    )

    public_vcf = subparsers.add_parser(
        "prepare-public-vcf",
        help="Validate a donor-matched public VCF and create phased per-sample seq/tag VCFs",
    )
    public_vcf.add_argument("--sample", required=True, help="Pipeline output sample ID")
    public_vcf.add_argument("--source-sample", required=True, help="Sample ID stored in the public VCF")
    public_vcf.add_argument("--source-vcf", required=True, action="append", type=Path, help="Repeat for chromosome- or type-split VCFs")
    public_vcf.add_argument("--reference-fasta", required=True, type=Path)
    public_vcf.add_argument("--output-dir", required=True, type=Path)
    public_vcf.add_argument("--callset", required=True)
    public_vcf.add_argument("--phase-source", choices=["provided"], default="provided")
    public_vcf.add_argument("--chroms", default="chr1-chr22")
    public_vcf.add_argument("--sequence-variant-policy", choices=["snp_only", "snv_indel"], default="snv_indel")
    public_vcf.add_argument("--bcftools-bin", default="bcftools")
    public_vcf.add_argument("--overwrite", action="store_true")

    personal_fasta = subparsers.add_parser(
        "create-personal-fastas",
        help="Create process_code-style HAP1/HAP2 personal FASTA and optional chain files from per-sample seq VCFs",
    )
    personal_fasta.add_argument("--cohort", required=True, help="Cohort/project label used in metadata")
    personal_fasta.add_argument("--seq-vcf-dir", required=True, type=Path, help="Directory containing <sample>.seq.vcf.gz from step 02")
    personal_fasta.add_argument("--output-root", required=True, type=Path, help="Final cohort root. Outputs go to <output-root>/<sample>/individual_fasta")
    personal_fasta.add_argument("--reference-fasta", required=True, type=Path, help="Personal-sequence consensus reference, e.g. hg38_autosomes.fa")
    personal_fasta.add_argument("--sample-list", type=Path, help="Text file with one sample ID per line")
    personal_fasta.add_argument("--sample", action="append", default=[], help="Sample ID. May be repeated; used when --sample-list is not provided")
    personal_fasta.add_argument("--chroms", default="chr1-chr22", help="Chromosomes to process, e.g. chr22 or chr1-chr22")
    personal_fasta.add_argument("--sequence-variant-policy", required=True, choices=["snp_only", "snv_indel"], help="Sequence FASTA policy")
    personal_fasta.add_argument("--bcftools-bin", default="bcftools")
    personal_fasta.add_argument("--samtools-bin", default="samtools")
    personal_fasta.add_argument("--jobs", type=int, default=8, help="Sample-level parallel jobs")
    personal_fasta.add_argument("--summary-tsv", type=Path)
    personal_fasta.add_argument("--metadata-json", type=Path)
    personal_fasta.add_argument("--sanitize-invalid-bases", action="store_true", help="Convert non-ACGTN FASTA bases to N instead of failing")
    personal_fasta.add_argument("--overwrite", action="store_true")

    haplotag = subparsers.add_parser(
        "create-haplotagged-alignments",
        help="Create process_code-style WhatsHap haplotagged BAMs and taglists from per-sample tag VCFs",
    )
    haplotag.add_argument("--cohort", required=True, help="Cohort/project label used in metadata")
    haplotag.add_argument("--raw-root", required=True, type=Path, help="Raw/source cohort root containing sample WGS alignments")
    haplotag.add_argument("--alignment", type=Path, help="Explicit alignment path for a single --sample; takes priority over --raw-root/--alignment-pattern")
    haplotag.add_argument("--tag-vcf-dir", required=True, type=Path, help="Directory containing <sample>.tag.vcf.gz from step 02")
    haplotag.add_argument("--output-root", required=True, type=Path, help="Final cohort root. Outputs go to <output-root>/<sample>/haplotagged")
    haplotag.add_argument("--reference-fasta", required=True, type=Path, help="CRAM decode/WhatsHap reference, e.g. Vienna")
    haplotag.add_argument("--platform", required=True, help="Sequencing platform; long-read public VCF requires explicit opt-in")
    haplotag.add_argument("--source-sample-id", help="Manifest의 명시적 donor alias; 단일 sample 실행에서만 사용")
    haplotag.add_argument(
        "--allow-longread-vcf-haplotag",
        action="store_true",
        help="Explicitly allow a prevalidated public phased VCF for long-read haplotagging",
    )
    haplotag.add_argument("--sample-list", type=Path, help="Text file with one sample ID per line")
    haplotag.add_argument("--sample", action="append", default=[], help="Sample ID. May be repeated; used when --sample-list is not provided")
    haplotag.add_argument("--alignment-pattern", default="{sample}/wgs/{sample}.hg38.cram", help="Path pattern under --raw-root. Supports {sample}")
    haplotag.add_argument("--output-prefix", default="{sample}.ont.hg38.haplotagged", help="Output BAM prefix pattern. Supports {sample}")
    haplotag.add_argument("--whatshap-bin", default="whatshap")
    haplotag.add_argument("--samtools-bin", default="samtools")
    haplotag.add_argument("--jobs", type=int, default=6, help="Sample-level parallel jobs")
    haplotag.add_argument("--output-threads", type=int, default=4, help="WhatsHap output/index threads per sample")
    haplotag.add_argument("--inspect-reads", type=int, default=50000, help="Number of reads to inspect for HP/PS before and after haplotagging")
    haplotag.add_argument("--tag-supplementary", action="store_true", help="Pass --tag-supplementary=copy-primary to WhatsHap")
    haplotag.add_argument("--summary-tsv", type=Path)
    haplotag.add_argument("--metadata-json", type=Path)
    haplotag.add_argument("--overwrite", action="store_true")

    clair3 = subparsers.add_parser(
        "call-longread-clair3",
        help="Call long-read variants with Clair3 and create HP/PS-tagged BAM with WhatsHap",
    )
    clair3.add_argument("--cohort", required=True, help="Cohort/project label used in metadata")
    clair3.add_argument("--output-root", required=True, type=Path, help="Final cohort root. Outputs go to <output-root>/<sample>/longread_calling and haplotagged")
    clair3.add_argument("--reference-fasta", required=True, type=Path, help="Reference used for CRAM decoding, Clair3 calling, and WhatsHap haplotagging")
    clair3.add_argument("--model-path", required=True, type=Path, help="Clair3 platform model directory")
    clair3.add_argument("--platform", required=True, help="Long-read platform. Accepted aliases resolve to Clair3 ont or hifi")
    clair3.add_argument("--raw-root", type=Path, help="Raw/source cohort root containing sample alignments")
    clair3.add_argument("--alignment", type=Path, help="Single-sample BAM/CRAM override. Requires exactly one --sample")
    clair3.add_argument("--sample-list", type=Path, help="Text file with one sample ID per line")
    clair3.add_argument("--sample", action="append", default=[], help="Sample ID. May be repeated; exactly one when --alignment is used")
    clair3.add_argument("--alignment-pattern", default="{sample}/wgs/{sample}.hg38.cram", help="Path pattern under --raw-root. Supports {sample}")
    clair3.add_argument("--output-prefix", default="{sample}.clair3_longread", help="Output prefix pattern. Supports {sample}")
    clair3.add_argument("--clair3-bin", default="run_clair3.py", help="Clair3 run_clair3.py path; .py entrypoints are executed with --python-bin")
    clair3.add_argument("--whatshap-bin", default="whatshap")
    clair3.add_argument("--samtools-bin", default="samtools")
    clair3.add_argument("--python-bin", default="python")
    clair3.add_argument("--pypy-bin", default="pypy3")
    clair3.add_argument("--parallel-bin", default="parallel")
    clair3.add_argument("--tabix-bin", default="tabix")
    clair3.add_argument("--jobs", type=int, default=1, help="Sample-level parallel jobs")
    clair3.add_argument("--threads", type=int, default=32, help="Clair3 threads per sample")
    clair3.add_argument("--output-threads", type=int, default=8, help="WhatsHap output/index threads per sample")
    clair3.add_argument("--inspect-reads", type=int, default=50000, help="Reads sampled for HP/PS validation after haplotagging")
    clair3.add_argument("--contigs", default="", help="Contigs passed to Clair3 ctg_name. Supports chr22 or chr1-chr22. Use this or --bed-fn")
    clair3.add_argument("--bed-fn", type=Path, help="Optional Clair3 BED restriction. Mutually exclusive with --contigs")
    clair3.add_argument("--chunk-size", type=int, default=5000000)
    clair3.add_argument("--min-mq", type=int, default=5)
    clair3.add_argument("--min-coverage", type=int, default=2)
    clair3.add_argument("--qual", type=int, default=2)
    clair3.add_argument("--enable-long-indel", action="store_true", help="Pass --enable_long_indel=True to Clair3")
    clair3.add_argument("--keep-intermediate-dir", action="store_true", help="Do not pass Clair3 --remove_intermediate_dir=True")
    clair3.add_argument("--no-tag-supplementary", action="store_true", help="Do not copy primary HP/PS to supplementary alignments in external WhatsHap haplotagging")
    clair3.add_argument("--ignore-read-groups", action="store_true", help="Pass --ignore-read-groups to WhatsHap haplotag")
    clair3.add_argument("--haplotag-mode", choices=["external_whatshap", "clair3_internal", "none"], default="external_whatshap", help="Default runs explicit WhatsHap haplotag after Clair3 phased VCF")
    clair3.add_argument("--export-sequence-vcf", action="store_true", help="Export the selected Clair3 VCF to <sequence-vcf-dir>/<sample>.seq.vcf.gz for long-read FASTA consensus")
    clair3.add_argument("--sequence-vcf-dir", type=Path, help="Directory where long-read-derived <sample>.seq.vcf.gz is written")
    clair3.add_argument("--sequence-vcf-export-mode", choices=["copy", "symlink"], default="copy", help="How to export the Clair3 VCF into the per-sample FASTA VCF slot")
    clair3.add_argument("--allow-unphased-sequence-vcf", action="store_true", help="Diagnostic only: allow FASTA consensus from a sequence VCF without phased heterozygous GTs")
    clair3.add_argument("--allow-empty-haplotag", action="store_true", help="Do not fail when haplotagged BAM has no HP tags in inspected reads")
    clair3.add_argument("--allow-pileup-vcf-fallback", action="store_true", help="Smoke-only: use pileup.vcf when Clair3 does not produce merge_output.vcf.gz")
    clair3.add_argument("--extra-clair3-arg", action="append", default=[], help="Additional raw Clair3 argument, e.g. --extra-clair3-arg=--pileup_only=True")
    clair3.add_argument("--summary-tsv", type=Path)
    clair3.add_argument("--metadata-json", type=Path)
    clair3.add_argument("--overwrite", action="store_true")
    clair3.add_argument("--dry-run", action="store_true", help="Write planned commands/metadata without executing Clair3 or WhatsHap")

    clair3_qc = subparsers.add_parser(
        "qc-longread-clair3",
        help="Inspect Clair3 VCF and HP/PS-tagged BAM quality for long-read samples",
    )
    clair3_qc.add_argument("--cohort", required=True)
    clair3_qc.add_argument("--output-root", required=True, type=Path)
    clair3_qc.add_argument("--platform", required=True, help="Long-read platform: ont or hifi")
    clair3_qc.add_argument("--model-path", required=True, type=Path, help="Selected Clair3 model directory")
    clair3_qc.add_argument("--reference-fasta", required=True, type=Path)
    clair3_qc.add_argument("--sample-list", type=Path)
    clair3_qc.add_argument("--sample", action="append", default=[])
    clair3_qc.add_argument("--samtools-bin", default="samtools")
    clair3_qc.add_argument("--inspect-reads", type=int, default=50000)
    clair3_qc.add_argument("--min-mq", type=int, default=5)
    clair3_qc.add_argument("--min-coverage", type=int, default=2)
    clair3_qc.add_argument("--qual", type=float, default=2.0)
    clair3_qc.add_argument("--min-hp-fraction", type=float, default=0.50)
    clair3_qc.add_argument("--output-prefix", default="{sample}.clair3_longread")
    clair3_qc.add_argument("--summary-tsv", type=Path)
    clair3_qc.add_argument("--metadata-json", type=Path)



    denovo = subparsers.add_parser(
        "run-denovo-assembly",
        help="Extract long reads from BAM/CRAM and run hifiasm --ont/HiFi de novo assembly benchmark",
    )
    denovo.add_argument("--cohort", required=True)
    denovo.add_argument("--output-root", required=True, type=Path)
    denovo.add_argument("--reference-fasta", required=True, type=Path, help="Reference for CRAM decoding during FASTQ extraction")
    denovo.add_argument("--platform", required=True, help="Long-read platform: ont or hifi")
    denovo.add_argument("--sample-list", type=Path)
    denovo.add_argument("--sample", action="append", default=[])
    denovo.add_argument("--raw-root", type=Path)
    denovo.add_argument("--alignment", type=Path, help="Single-sample BAM/CRAM override. Requires exactly one --sample")
    denovo.add_argument("--alignment-pattern", default="{sample}/wgs/{sample}.hg38.cram")
    denovo.add_argument("--contigs", default="chr1-chr22", help="Contigs to extract for benchmark/assembly, e.g. chr22 or chr1-chr22")
    denovo.add_argument("--label", default="", help="Output label under <sample>/denovo_assembly. Default is assembler_platform_scope")
    denovo.add_argument("--assembler", choices=["hifiasm"], default="hifiasm")
    denovo.add_argument("--hifiasm-bin", default="hifiasm")
    denovo.add_argument("--samtools-bin", default="samtools")
    denovo.add_argument("--threads", type=int, default=48)
    denovo.add_argument("--fastq-threads", type=int, default=8)
    denovo.add_argument("--exclude-flags", type=int, default=2308, help="Default excludes unmapped, secondary, supplementary alignments")
    denovo.add_argument("--min-contig-len", type=int, default=100000)
    denovo.add_argument("--keep-fastq", dest="keep_fastq", action="store_true", default=True)
    denovo.add_argument("--no-keep-fastq", dest="keep_fastq", action="store_false")
    denovo.add_argument("--keep-gfa", dest="keep_gfa", action="store_true", default=True)
    denovo.add_argument("--no-keep-gfa", dest="keep_gfa", action="store_false")
    denovo.add_argument("--keep-bin", dest="keep_bin", action="store_true", default=False, help="Keep hifiasm *.bin restart/intermediate files")
    denovo.add_argument("--no-keep-bin", dest="keep_bin", action="store_false")
    denovo.add_argument("--overwrite", action="store_true")
    denovo.add_argument("--dry-run", action="store_true")
    denovo.add_argument("--summary-tsv", type=Path)
    denovo.add_argument("--metadata-json", type=Path)

    anchor = subparsers.add_parser(
        "anchor-denovo-assembly",
        help="Align de novo assembly FASTA contigs to a reference and write anchoring QC",
    )
    anchor.add_argument("--cohort", required=True)
    anchor.add_argument("--output-root", required=True, type=Path)
    anchor.add_argument("--reference-fasta", required=True, type=Path)
    anchor.add_argument("--sample", required=True)
    anchor.add_argument("--assembly-dir", type=Path)
    anchor.add_argument("--label", default="")
    anchor.add_argument("--minimap2-bin", default="minimap2")
    anchor.add_argument("--preset", default="asm20")
    anchor.add_argument("--threads", type=int, default=48)
    anchor.add_argument("--min-mapq", type=int, default=30)
    anchor.add_argument("--min-contig-len", type=int, default=100000)
    anchor.add_argument("--overwrite", action="store_true")
    anchor.add_argument("--dry-run", action="store_true")
    anchor.add_argument("--summary-tsv", type=Path)
    anchor.add_argument("--metadata-json", type=Path)

    bamonly_phase = subparsers.add_parser(
        "phase-longread-bamonly",
        help="Infer unoriented long-read BAM-only read clusters/phase blocks without VCF haplotagging",
    )
    bamonly_phase.add_argument("--sample-id", required=True)
    bamonly_phase.add_argument("--alignment", required=True, type=Path)
    bamonly_phase.add_argument("--output-dir", required=True, type=Path)
    bamonly_phase.add_argument("--reference-fasta", type=Path, help="Reference FASTA for CRAM decoding/pileup, e.g. Vienna")
    bamonly_phase.add_argument("--contigs", default="chr1-chr22", help="Contigs to scan when --regions is not given")
    bamonly_phase.add_argument("--regions", default="", help="Optional 1-based regions, e.g. chr22:20000000-21000000. Overrides whole-contig scanning")
    bamonly_phase.add_argument("--platform", default="ont")
    bamonly_phase.add_argument("--min-mapq", type=int, default=20)
    bamonly_phase.add_argument("--min-base-quality", type=int, default=0, help="Default 0 because some ONT CRAMs do not retain usable per-base quality for count_coverage filtering")
    bamonly_phase.add_argument("--min-depth", type=int, default=6, help="Minimum valid A/C/G/T read depth for a candidate marker")
    bamonly_phase.add_argument("--min-marker-reads", type=int, default=2, help="Minimum supporting reads for the minor allele")
    bamonly_phase.add_argument("--min-minor-allele-fraction", type=float, default=0.20)
    bamonly_phase.add_argument("--max-major-allele-fraction", type=float, default=0.80)
    bamonly_phase.add_argument("--max-third-allele-fraction", type=float, default=0.10)
    bamonly_phase.add_argument("--max-marker-gap", type=int, default=1_000_000, help="Split phase-block candidates when adjacent markers are farther apart")
    bamonly_phase.add_argument("--max-markers-per-block", type=int, default=5_000)
    bamonly_phase.add_argument("--min-read-markers", type=int, default=2, help="Minimum marker observations required to cluster a read")
    bamonly_phase.add_argument("--min-score-delta", type=int, default=2, help="Minimum score difference between cluster profiles")
    bamonly_phase.add_argument("--max-iterations", type=int, default=6)
    bamonly_phase.add_argument("--include-supplementary", action="store_true", help="Use supplementary alignments for marker discovery/clustering")
    bamonly_phase.add_argument("--exclude-duplicates", action="store_true", help="Exclude duplicate reads from marker discovery/clustering")
    bamonly_phase.add_argument("--max-pileup-depth", type=int, default=100_000, help="Kept for compatibility; count_coverage marker detection does not use pileup max-depth directly")
    bamonly_phase.add_argument("--marker-detection-chunk-bp", type=int, default=1_000_000, help="Chunk size for C-backed A/C/G/T coverage counting during marker detection")
    bamonly_phase.add_argument("--max-markers-per-contig", type=int, default=0, help="0 means unlimited; intended for smoke tests/debugging")
    bamonly_phase.add_argument("--unknown-read-policy", choices=["include_both", "drop_unknown"], default="include_both", help="Recorded for downstream routing; cluster=0 rows are still emitted in this step")
    bamonly_phase.add_argument("--dry-run", action="store_true", help="Print the planned contract without scanning the BAM/CRAM")

    personal_realign = subparsers.add_parser(
        "create-personal-realigned-alignments",
        help="Create process_code 06-style HAP1/HAP2 personal-coordinate BAM/CRAM alignments",
    )
    personal_realign.add_argument("--cohort", required=True, help="Cohort/project label used in summary files")
    personal_realign.add_argument("--output-root", required=True, type=Path, help="Final cohort root containing <sample>/individual_fasta")
    personal_realign.add_argument("--sample-list", type=Path, help="Text file with one sample ID per line")
    personal_realign.add_argument("--sample", action="append", default=[], help="Sample ID. May be repeated; required unless --sample-list is used")
    personal_realign.add_argument("--chroms", default="chr1-chr22", help="Chromosomes to realign, e.g. chr22 or chr1-chr22")
    personal_realign.add_argument("--sequence-variant-policy", required=True, choices=["snp_only", "snv_indel"], help="snv_indel is the normal pre-training use case")
    personal_realign.add_argument("--sequence-source", choices=["vcf_consensus", "assembly_direct"], default="vcf_consensus")
    personal_realign.add_argument("--assembly-phase-labels-validated", action="store_true", help="assembly_direct의 HP 기반 HAP routing에는 VCF와의 HAP orientation 검증이 필요함")
    personal_realign.add_argument("--dataset-role", required=True, choices=["pre_training", "post_training"], help="06 is pre_training-only by default")
    personal_realign.add_argument("--alignment", type=Path, help="Single-sample haplotagged hg38 BAM/CRAM override")
    personal_realign.add_argument("--alignment-pattern", default="{sample}/haplotagged/{sample}.ont.hg38.haplotagged.bam", help="Pattern relative to --output-root. Supports {sample}")
    personal_realign.add_argument("--hap1-fasta-dir", type=Path, help="Single-sample HAP1 FASTA directory override")
    personal_realign.add_argument("--hap2-fasta-dir", type=Path, help="Single-sample HAP2 FASTA directory override")
    personal_realign.add_argument("--output-prefix", default="{sample}.personal.HAP{hap}", help="Output prefix pattern. Supports {sample}, {hap}, {haplotype}")
    personal_realign.add_argument("--realigned-subdir", default="realigned", help="Relative directory below the sample root; use a HAP-mode-specific path to prevent output collisions")
    personal_realign.add_argument("--reference-fasta", type=Path, help="Reference FASTA for input CRAM decoding, if --alignment or pattern points to CRAM")
    personal_realign.add_argument("--format", dest="output_format", choices=["cram"], default="cram", help="CRAM only: Stage 12 (legacy Stage70) verifies target.fa against CRAM @SQ M5")
    personal_realign.add_argument("--platform", default="", help="Read platform used to resolve --preset auto: ont->map-ont, hifi/ccs->map-hifi, clr->map-pb, illumina->sr")
    personal_realign.add_argument("--preset", default="auto", help="minimap2 preset. Default auto resolves from --platform.")
    personal_realign.add_argument("--threads", type=int, default=24, help="aligner threads per haplotype")
    personal_realign.add_argument("--sort-threads", type=int, default=8, help="samtools sort threads per haplotype")
    personal_realign.add_argument("--sort-mem", default="2G", help="samtools sort memory per thread")
    personal_realign.add_argument("--jobs", type=int, default=1, help="Sample-level parallelism")
    personal_realign.add_argument("--hap-jobs", type=int, choices=[1, 2], default=1, help="HAP1/HAP2 parallelism within each sample")
    personal_realign.add_argument("--route-jobs", type=int, default=1, help="Indexed contig workers for Illumina all_to_both routing")
    personal_realign.add_argument("--min-mapped-frac", type=float, default=0.90, help="Legacy routed-mapped reference; enforced only with --ratio-qc-policy strict")
    personal_realign.add_argument(
        "--ratio-qc-policy", choices=("record_only", "strict"), default="record_only",
        help="record_only stores mapping/pair ratios without gating; strict reproduces legacy hard cutoffs",
    )
    personal_realign.add_argument("--hap-read-mode", choices=["all_to_both", "hp_routed_with_unknown", "hap_specific"], default="", help="High-level routing mode. all_to_both preserves coverage; hp_routed_with_unknown duplicates HP0/missing; hap_specific retains only HP-consistent reads/fragments. When set, this mode overrides legacy routing flags.")
    personal_realign.add_argument("--read-routing-policy", choices=["hp_routed", "all_to_both"], default="hp_routed", help="Legacy low-level routing; use --hap-read-mode in current Stage 10.")
    personal_realign.add_argument("--unknown-read-policy", choices=["include_both", "drop_unknown"], default="include_both", help="How hp_routed handles HP0/missing reads: include_both routes to both HAPs, drop_unknown excludes them. Default: include_both")
    personal_realign.add_argument("--include-unmapped", dest="include_unmapped", action="store_true", default=True, help="Include original unmapped reads with XU:i:1 tag")
    personal_realign.add_argument("--no-include-unmapped", dest="include_unmapped", action="store_false", help="Do not include original unmapped reads")
    personal_realign.add_argument("--verify-counts", dest="verify_counts", action="store_true", default=True, help="Verify routed/rescue primary counts against input counts")
    personal_realign.add_argument("--no-verify-counts", dest="verify_counts", action="store_false", help="Skip exact count verification")
    personal_realign.add_argument("--count-tolerance-abs", type=int, default=10, help="Maximum missing primary records allowed per count check; extra records still fail")
    personal_realign.add_argument("--keep-target", action="store_true", help="Compatibility option; CRAM output always keeps target.fa")
    personal_realign.add_argument("--allow-non-pretraining", action="store_true", help="Validation-only override for post_training or snp_only")
    personal_realign.add_argument("--samtools-bin", default="samtools")
    personal_realign.add_argument("--minimap2-bin", default="minimap2")
    personal_realign.add_argument("--bwa-mem2-bin", default="bwa-mem2")
    personal_realign.add_argument("--illumina-aligner", choices=["bwa-mem2", "minimap2"], default="bwa-mem2", help="Illumina personal-HAP aligner; BWA-MEM2 is the production default")
    personal_realign.add_argument("--bwa-mem2-index-root", type=Path, help="Reference-fingerprint BWA-MEM2 index cache")
    personal_realign.add_argument("--illumina-duplicate-mode", choices=["preserve_source", "none"], default="preserve_source", help="Transport source duplicate flag through FASTQ via OD tag, or clear it")
    personal_realign.add_argument("--inspect-reads", type=int, default=5000, help="Reads sampled for NM/base comparison")
    personal_realign.add_argument("--dry-run", action="store_true", help="Build target/scripts/metadata but do not run minimap2")
    personal_realign.add_argument("--overwrite", action="store_true")
    personal_realign.add_argument("--summary-tsv", type=Path)
    personal_realign.add_argument("--metadata-json", type=Path)

    window_index = subparsers.add_parser(
        "create-window-index",
        help="Create process_code-style HAP1/HAP2 window index files with integrated deterministic random offsets",
    )
    window_index.add_argument("--cohort", required=True, help="Cohort/project label used in metadata")
    window_index.add_argument("--output-root", required=True, type=Path, help="Final cohort root containing <sample>/individual_fasta")
    window_index.add_argument("--sample-list", type=Path, help="Text file with one sample ID per line")
    window_index.add_argument("--sample", action="append", default=[], help="Sample ID. May be repeated; used when --sample-list is not provided")
    window_index.add_argument("--chroms", default="chr1-chr22", help="Chromosomes to process, e.g. chr22 or chr1-chr22")
    window_index.add_argument("--sequence-variant-policy", required=True, choices=["snp_only", "snv_indel"], help="Personal FASTA policy; snv_indel requires chain files")
    window_index.add_argument("--sequence-source", choices=["vcf_consensus", "assembly_direct"], default="vcf_consensus")
    window_index.add_argument("--window-size", type=int, default=1 << 20)
    window_index.add_argument("--offset-unit", type=int, help="Deterministic random offset range. Default: --window-size")
    window_index.add_argument("--fixed-window-start", action="store_true", help="Use process_code-style 0-start windows by setting offset_unit=1")
    window_index.add_argument("--seed", type=int, default=42)
    window_index.add_argument("--drop-n", action=argparse.BooleanOptionalAction, default=False, help="Legacy override: drop any non-ACGT; default preserves internal N and excludes all-N")
    window_index.add_argument("--blacklist-bed", type=Path, default=Path(DEFAULT_BLACKLIST_BED), help="Reference-coordinate blacklist BED")
    window_index.add_argument("--disable-reference-blacklist", action="store_true", help="assembly_direct에서 GRCh38 multi-chain projection 전 reference blacklist를 비활성화")
    window_index.add_argument("--max-blacklist-frac", type=float, default=1.0, help="Default records blacklist overlap without removing windows")
    window_index.add_argument("--jobs", type=int, default=8, help="Sample-level parallel jobs")
    window_index.add_argument("--summary-tsv", type=Path)
    window_index.add_argument("--metadata-json", type=Path)
    window_index.add_argument("--overwrite", action="store_true")

    full = subparsers.add_parser(
        "preprocess-full-shared-hap",
        help="Run full WGS shared-HAP preprocessing, compact to strict 4bit A/C/G/T/N/P, and write metadata",
    )
    full.add_argument("--sample-id", required=True)
    full.add_argument("--source-sample-id", default="")
    full.add_argument("--dataset-role", required=True, choices=["pre_training", "post_training"])
    full.add_argument("--read-type", required=True, choices=["long_read", "short_read"])
    full.add_argument("--project-id", required=True)
    full.add_argument("--project-label", default="")
    full.add_argument("--platform", default="")
    full.add_argument("--reference-build", default="hg38")
    full.add_argument("--sample-status", default="healthy")
    full.add_argument("--disease", default="NaN")
    full.add_argument("--hap1-fasta-dir", type=Path)
    full.add_argument("--hap1-chain-dir", type=Path)
    full.add_argument("--hap2-fasta-dir", type=Path)
    full.add_argument("--hap2-chain-dir", type=Path)
    full.add_argument("--alignment", required=True, type=Path)
    full.add_argument("--contigs", required=True, help="Comma/range contigs: chr22,chr1-chr22,1-22")
    full.add_argument("--dataset-root", type=Path)
    full.add_argument("--output-sample-root", type=Path)
    full.add_argument("--work-root", type=Path)
    full.add_argument("--window-index-dir", type=Path, help="Optional <sample>/windows directory from create-window-index; kept windows are reused by direct compact")
    full.add_argument("--read-coordinate-mode", choices=["reference_projected", "personal_realigned"], default="personal_realigned", help="personal_realigned is the v2 default and uses HAP-specific realigned alignments/split readsets; reference_projected keeps the legacy shared readset")
    full.add_argument("--hap1-alignment", type=Path, help="HAP1 personal-coordinate BAM/CRAM for --read-coordinate-mode personal_realigned")
    full.add_argument("--hap2-alignment", type=Path, help="HAP2 personal-coordinate BAM/CRAM for --read-coordinate-mode personal_realigned")
    full.add_argument("--layout-name", default=FULL_LAYOUT_NAME, help="Compact layout directory name; mode-specific names prevent incompatible HAP-routing outputs from colliding")
    full.add_argument("--hap-read-mode", choices=["all_to_both", "hp_routed_with_unknown", "hap_specific"], default="all_to_both")
    full.add_argument("--storage-policy", choices=["model_input_only", "keep_raw_intermediate"], default="model_input_only")
    full.add_argument("--jobs", type=int, default=1)
    full.add_argument("--windows-per-chunk", type=int, default=32)
    full.add_argument("--compact-mode", choices=["direct", "two_step"], default="direct", help="direct writes final 4bit compact zarr without raw object zarr; two_step keeps the legacy raw zarr conversion path")
    full.add_argument("--validation-level", choices=["full", "standard", "fast"], default="full")
    full.add_argument("--compact-compressor", choices=sorted(SUPPORTED_COMPACT_COMPRESSORS), default=DEFAULT_COMPACT_COMPRESSOR, help="zarr chunk compressor for direct compact output")
    full.add_argument("--compact-compression-level", type=int, default=DEFAULT_COMPACT_COMPRESSION_LEVEL, help="Blosc compression level 0..9 for direct compact output")
    full.add_argument("--profile-direct-compact", action="store_true", help="write detailed per-chunk direct compact timing; intended for profiling only")
    full.add_argument("--schedule", choices=["input", "longest-first"], default="longest-first")
    full.add_argument("--window-size", type=int, default=1 << 20)
    full.add_argument("--offset-unit", type=int)
    full.add_argument("--fixed-window-start", action="store_true", help="Use fixed 0-start windows by setting offset_unit=1; default uses deterministic random offset within the full window size")
    full.add_argument("--seed", type=int, default=42)
    full.add_argument("--min-mapq", type=int, default=30)
    full.add_argument("--stage70-filter-preset", choices=["auto", "short_read_base_v1", "long_read_primary_supp_evidence_v1"], default="auto", help="auto resolves short/long Stage 12 compact filter (legacy Stage70 schema) and primary/supplementary depth policy")
    full.add_argument("--max-windows", type=int)
    full.add_argument("--overwrite", action="store_true")
    full.add_argument("--no-resume", action="store_true")
    full.add_argument("--allow-untagged", action="store_true")
    full.add_argument("--haplotag-if-missing", action="store_true")
    full.add_argument("--allow-unphased-longread", action="store_true", help="Ablation/debug override allowing long-read compact export with no HP tags. Production default blocks it.")
    full.add_argument("--auto-prepare-hap-fasta", action="store_true")
    full.add_argument("--phased-vcf", type=Path)
    full.add_argument("--common-vcf-dir", type=Path)
    full.add_argument("--common-vcf-pattern", default="")
    full.add_argument("--common-vcf-filter-expr", default='GT="het"')
    full.add_argument("--consensus-filter-expr", default='(FILTER="PASS" || FILTER=".")')
    full.add_argument("--haplotag-vcf-types", default="snps,indels")
    full.add_argument("--sequence-variant-policy", choices=["snp_only", "snv_indel"], default="snv_indel")
    full.add_argument("--reference-fasta", type=Path, help="Reference used for CRAM decoding and WhatsHap haplotagging, e.g. Vienna")
    full.add_argument("--sequence-reference-fasta", type=Path, help="Reference used for bcftools consensus personal FASTA, e.g. UCSC hg38 autosomes")
    full.add_argument("--personal-fasta-qc", type=Path, help="Stage 03 personal_fasta_qc.json; personal_realigned mode validates exact reference/chain provenance")
    full.add_argument("--phase-inspect-max-reads", type=int, default=50000)
    full.add_argument("--synthesize-missing-ps", action="store_true", help="Optionally convert missing/zero PS on HP-tagged reads to synthetic PS=1; default preserves missing/zero PS")
    full.add_argument("--whatshap-bin", default="whatshap")
    full.add_argument("--samtools-bin", default="samtools")
    full.add_argument("--bcftools-bin", default="bcftools")
    full.add_argument("--source-wgs", action="append", type=Path, default=[])
    full.add_argument("--source-vcf", action="append", type=Path, default=[])
    full.add_argument("--rna-bigwig", action="append", type=Path, default=[])

    igsr_manifest = subparsers.add_parser(
        "build-igsr-ont-manifest",
        help="Build preprocess-full-manifest TSV for IGSR ONT WGS-only CRAM directories",
    )
    igsr_manifest.add_argument("--raw-root", required=True, type=Path)
    igsr_manifest.add_argument("--output-manifest", required=True, type=Path)
    igsr_manifest.add_argument("--dataset-root", required=True, type=Path)
    igsr_manifest.add_argument("--common-vcf-dir", required=True, type=Path)
    igsr_manifest.add_argument("--reference-fasta", required=True, type=Path, help="Reference used for CRAM decoding and WhatsHap haplotagging, e.g. Vienna")
    igsr_manifest.add_argument("--contigs", default="chr1-chr22")
    igsr_manifest.add_argument("--dataset-role", default="pre_training", choices=["pre_training", "post_training"])
    igsr_manifest.add_argument("--read-type", default="long_read", choices=["long_read", "short_read"])
    igsr_manifest.add_argument("--project-id", default="igsr_ont_wgs_only")
    igsr_manifest.add_argument("--project-label", default="IGSR ONT WGS-only")
    igsr_manifest.add_argument("--platform", default="OxfordNanopore")
    igsr_manifest.add_argument("--reference-build", default="hg38")
    igsr_manifest.add_argument("--sample-status", default="healthy")
    igsr_manifest.add_argument("--disease", default="NaN")
    igsr_manifest.add_argument("--storage-policy", choices=["model_input_only", "keep_raw_intermediate"], default="model_input_only")
    igsr_manifest.add_argument("--jobs", type=int, default=32)
    igsr_manifest.add_argument("--windows-per-chunk", type=int, default=32)
    igsr_manifest.add_argument("--compact-mode", choices=["direct", "two_step"], default="direct")
    igsr_manifest.add_argument("--read-coordinate-mode", choices=["reference_projected", "personal_realigned"], default="personal_realigned", help="Manifest-level read coordinate contract for final compact. v2 default: personal_realigned")
    igsr_manifest.add_argument("--hap1-alignment-pattern", default="{sample}/realigned/HAP1/{sample}.personal.HAP1.cram", help="Pattern under dataset_root/final/<role>/<read_type>/<project_id>. Supports {sample}, {hap}, {haplotype}")
    igsr_manifest.add_argument("--hap2-alignment-pattern", default="{sample}/realigned/HAP2/{sample}.personal.HAP2.cram", help="Pattern under dataset_root/final/<role>/<read_type>/<project_id>. Supports {sample}, {hap}, {haplotype}")
    igsr_manifest.add_argument("--haplotagged-alignment-pattern", default="{sample}/haplotagged/{sample}.clair3_longread.haplotagged.bam", help="Haplotagged BAM pattern used as phase-inspection/source alignment in personal_realigned final compact mode. v2 long-read default is Clair3 + WhatsHap")
    igsr_manifest.add_argument("--window-index-pattern", default="{sample}/windows", help="Window index directory pattern used in personal_realigned final compact mode")
    igsr_manifest.add_argument("--validation-level", choices=["full", "standard", "fast"], default="full")
    igsr_manifest.add_argument("--compact-compressor", choices=sorted(SUPPORTED_COMPACT_COMPRESSORS), default=DEFAULT_COMPACT_COMPRESSOR, help="zarr chunk compressor for direct compact output")
    igsr_manifest.add_argument("--compact-compression-level", type=int, default=DEFAULT_COMPACT_COMPRESSION_LEVEL, help="Blosc compression level 0..9 for direct compact output")
    igsr_manifest.add_argument("--schedule", choices=["input", "longest-first"], default="longest-first")
    igsr_manifest.add_argument("--window-size", type=int, default=1 << 20)
    igsr_manifest.add_argument("--offset-unit", type=int)
    igsr_manifest.add_argument("--fixed-window-start", action="store_true", help="Use fixed 0-start windows by setting offset_unit=1; default uses deterministic random offset within the full window size")
    igsr_manifest.add_argument("--seed", type=int, default=42)
    igsr_manifest.add_argument("--min-mapq", type=int, default=30)
    igsr_manifest.add_argument("--max-windows", type=int)
    igsr_manifest.add_argument("--phase-inspect-max-reads", type=int, default=50000)
    igsr_manifest.add_argument("--whatshap-bin", default="whatshap")
    igsr_manifest.add_argument("--samtools-bin", default="samtools")
    igsr_manifest.add_argument("--bcftools-bin", default="bcftools")
    igsr_manifest.add_argument("--common-vcf-pattern", default="")
    igsr_manifest.add_argument("--common-vcf-filter-expr", default='GT="het"')
    igsr_manifest.add_argument("--consensus-filter-expr", default='(FILTER="PASS" || FILTER=".")')
    igsr_manifest.add_argument("--haplotag-vcf-types", default="snps,indels")
    igsr_manifest.add_argument("--sequence-variant-policy", choices=["snp_only", "snv_indel"], default="snv_indel")
    igsr_manifest.add_argument("--sequence-reference-fasta", required=True, type=Path, help="UCSC hg38 reference used for personal HAP FASTA consensus")
    igsr_manifest.add_argument("--synthesize-missing-ps", action="store_true")
    igsr_manifest.add_argument("--allow-unphased-longread", action="store_true", help="Write manifest rows that allow unphased long-read ablation runs")
    igsr_manifest.add_argument("--sample", action="append", default=[])
    igsr_manifest.add_argument("--sample-file", type=Path)
    igsr_manifest.add_argument("--pilot-size", type=int)
    igsr_manifest.add_argument("--limit", type=int)
    igsr_manifest.add_argument("--summary-tsv", type=Path)
    igsr_manifest.add_argument("--overwrite-manifest", action="store_true")

    validate = subparsers.add_parser(
        "validate-compact-v2",
        help="Validate v2 compact contig layout, including hap_sequence/readset/depth consistency",
    )
    validate.add_argument("--compact-contig-dir", required=True, type=Path)
    validate.add_argument("--read-coordinate-mode", choices=["personal_realigned", "reference_projected"], default="personal_realigned")
    validate.add_argument("--validation-level", choices=["fast", "standard", "full"], default="standard")
    validate.add_argument("--output-json", type=Path)

    trace = subparsers.add_parser(
        "validate-stage70-trace",
        help="Independently compare personal HAP source CRAM records against existing Stage 12 compact Zarr (legacy Stage70 schema)",
    )
    trace.add_argument("--compact-root", required=True, type=Path, help="Existing compact directory containing chr*/")
    trace.add_argument("--contigs", required=True, help="Comma-separated contigs to validate")
    trace.add_argument("--report-root", required=True, type=Path, help="Separate QC directory; cannot be inside compact-root")
    trace.add_argument("--workers", type=int, default=1, help="Flat contig×HAP process count")
    trace.add_argument("--max-examples", type=int, default=20)
    trace.add_argument("--scratch-dir", type=Path, help="Optional local scratch for temporary SQLite multisets")
    full_manifest = subparsers.add_parser(
        "preprocess-full-manifest",
        help="Run full shared-HAP preprocessing for sample rows in a TSV manifest",
    )
    full_manifest.add_argument("--manifest", required=True, type=Path)
    full_manifest.add_argument("--status-tsv", required=True, type=Path)
    full_manifest.add_argument("--dataset-root", type=Path)
    full_manifest.add_argument("--sample-jobs", type=int, default=1)
    full_manifest.add_argument("--jobs", type=int, default=1, help="contig jobs per sample")
    full_manifest.add_argument("--overwrite", action="store_true")
    full_manifest.add_argument("--no-resume", action="store_true")
    full_manifest.add_argument("--storage-policy", choices=["model_input_only", "keep_raw_intermediate"], default="model_input_only")
    full_manifest.add_argument("--sample-catalog", type=Path, help="Optional sample_catalog.tsv to upsert after each completed sample")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "validate-wgs-input-manifest":
        contract = validate_wgs_input_manifest(
            args.manifest,
            args.sample,
            read_files_manifest=args.read_files_manifest,
            allow_missing_files=args.allow_missing_files,
        )
        if args.output_json:
            write_contract_json(contract, args.output_json)
        print(json.dumps(contract, ensure_ascii=False, indent=2))
        return
    if args.command == "prototype-readset":
        config = PrototypeConfig(
            sample_id=args.sample_id,
            alignment=args.alignment,
            contig=args.contig,
            output_dir=args.output_dir,
            min_mapq=args.min_mapq,
            max_reads=args.max_reads,
            exclude_duplicate=not args.keep_duplicates,
            exclude_supplementary=args.exclude_supplementary,
        )
        run_prototype(config)
        return
    if args.command == "convert-legacy-shared-hap":
        contigs = None
        if args.contigs:
            contigs = [item.strip() for item in args.contigs.split(",") if item.strip()]
        config = LegacyConvertConfig(
            sample_id=args.sample_id,
            source_sample_id=args.source_sample_id,
            legacy_root=args.legacy_root,
            output_sample_root=args.output_sample_root,
            dataset_role=args.dataset_role,
            read_type=args.read_type,
            project_id=args.project_id,
            contigs=contigs,
            layout_name=args.layout_name,
            overwrite=args.overwrite,
            chunk_records=args.chunk_records,
            hap_invalid_policy=args.hap_invalid_policy,
            rna_bigwigs=args.rna_bigwig,
            source_wgs=args.source_wgs,
            source_vcfs=args.source_vcf,
        )
        run_legacy_convert(config)
        return

    if args.command == "convert-legacy-manifest":
        config = ManifestRunConfig(
            manifest=args.manifest,
            status_tsv=args.status_tsv,
            overwrite=args.overwrite,
            jobs=args.jobs,
            hap_jobs=args.hap_jobs,
        )
        run_manifest(config)
        return

    if args.command == "convert-legacy-manifest-contigs":
        config = ManifestContigRunConfig(
            manifest=args.manifest,
            status_tsv=args.status_tsv,
            overwrite=args.overwrite,
            jobs=args.jobs,
        )
        run_manifest_contigs(config)
        return

    if args.command == "extract-cohort-common-vcf":
        config = CohortVcfConfig(
            cohort=args.cohort,
            sample_list=args.sample_list,
            samples=args.sample,
            output_dir=args.output_dir,
            chroms=[args.chroms],
            panel_base=args.panel_base,
            panel_file_template=args.panel_file_template,
            bcftools_bin=args.bcftools_bin,
            filter_status=args.filter_status,
            min_ac=args.min_ac,
            jobs=args.jobs,
            threads_per_chrom=args.threads_per_chrom,
            overwrite=args.overwrite,
            summary_tsv=args.summary_tsv,
            metadata_json=args.metadata_json,
            idx_cache_dir=args.idx_cache_dir,
        )
        summary = run_cohort_vcf_extract(config)
        print(summary)
        return

    if args.command == "extract-per-sample-vcfs":
        config = PerSampleVcfConfig(
            extraction_mode=args.extraction_mode,
            cohort=args.cohort,
            common_vcf_dir=args.common_vcf_dir,
            output_dir=args.output_dir,
            sample_list=args.sample_list,
            samples=args.sample,
            chroms=[args.chroms],
            sequence_variant_policy=args.sequence_variant_policy,
            common_vcf_pattern=args.common_vcf_pattern,
            concat_dir=args.concat_dir,
            bcftools_bin=args.bcftools_bin,
            jobs=args.jobs,
            concat_threads=args.concat_threads,
            overwrite=args.overwrite,
            fast_validate=args.fast_validate,
            summary_tsv=args.summary_tsv,
            metadata_json=args.metadata_json,
        )
        summary = run_per_sample_vcf_extract(config)
        print(summary)
        return

    if args.command == "prepare-public-vcf":
        config = PublicVcfConfig(
            sample=args.sample,
            source_sample=args.source_sample,
            source_vcfs=args.source_vcf,
            reference_fasta=args.reference_fasta,
            output_dir=args.output_dir,
            callset=args.callset,
            phase_source=args.phase_source,
            chroms=[args.chroms],
            sequence_variant_policy=args.sequence_variant_policy,
            bcftools_bin=args.bcftools_bin,
            overwrite=args.overwrite,
        )
        summary = prepare_public_vcf(config)
        print(summary)
        return

    if args.command == "create-personal-fastas":
        config = PersonalFastaConfig(
            cohort=args.cohort,
            seq_vcf_dir=args.seq_vcf_dir,
            output_root=args.output_root,
            reference_fasta=args.reference_fasta,
            sample_list=args.sample_list,
            samples=args.sample,
            chroms=[args.chroms],
            sequence_variant_policy=args.sequence_variant_policy,
            bcftools_bin=args.bcftools_bin,
            samtools_bin=args.samtools_bin,
            jobs=args.jobs,
            overwrite=args.overwrite,
            sanitize_invalid_bases=args.sanitize_invalid_bases,
            summary_tsv=args.summary_tsv,
            metadata_json=args.metadata_json,
        )
        summary = run_personal_fasta_generation(config)
        print(summary)
        return

    if args.command == "create-haplotagged-alignments":
        config = HaplotagConfig(
            cohort=args.cohort,
            raw_root=args.raw_root,
            alignment=args.alignment,
            tag_vcf_dir=args.tag_vcf_dir,
            output_root=args.output_root,
            reference_fasta=args.reference_fasta,
            platform=args.platform,
            allow_longread_vcf_haplotag=args.allow_longread_vcf_haplotag,
            source_sample_id=args.source_sample_id,
            sample_list=args.sample_list,
            samples=args.sample,
            alignment_pattern=args.alignment_pattern,
            output_prefix=args.output_prefix,
            whatshap_bin=args.whatshap_bin,
            samtools_bin=args.samtools_bin,
            jobs=args.jobs,
            output_threads=args.output_threads,
            inspect_reads=args.inspect_reads,
            tag_supplementary=args.tag_supplementary,
            overwrite=args.overwrite,
            summary_tsv=args.summary_tsv,
            metadata_json=args.metadata_json,
        )
        summary = run_haplotag_generation(config)
        print(summary)
        return


    if args.command == "call-longread-clair3":
        config = Clair3LongReadConfig(
            cohort=args.cohort,
            output_root=args.output_root,
            reference_fasta=args.reference_fasta,
            model_path=args.model_path,
            platform=args.platform,
            raw_root=args.raw_root,
            alignment=args.alignment,
            sample_list=args.sample_list,
            samples=args.sample,
            alignment_pattern=args.alignment_pattern,
            output_prefix=args.output_prefix,
            clair3_bin=args.clair3_bin,
            whatshap_bin=args.whatshap_bin,
            samtools_bin=args.samtools_bin,
            python_bin=args.python_bin,
            pypy_bin=args.pypy_bin,
            parallel_bin=args.parallel_bin,
            tabix_bin=args.tabix_bin,
            jobs=args.jobs,
            threads=args.threads,
            output_threads=args.output_threads,
            inspect_reads=args.inspect_reads,
            contigs=args.contigs,
            bed_fn=args.bed_fn,
            chunk_size=args.chunk_size,
            min_mq=args.min_mq,
            min_coverage=args.min_coverage,
            qual=args.qual,
            enable_long_indel=args.enable_long_indel,
            remove_intermediate_dir=not args.keep_intermediate_dir,
            tag_supplementary=not args.no_tag_supplementary,
            ignore_read_groups=args.ignore_read_groups,
            haplotag_mode=args.haplotag_mode,
            export_sequence_vcf=args.export_sequence_vcf,
            sequence_vcf_dir=args.sequence_vcf_dir,
            sequence_vcf_export_mode=args.sequence_vcf_export_mode,
            allow_unphased_sequence_vcf=args.allow_unphased_sequence_vcf,
            allow_empty_haplotag=args.allow_empty_haplotag,
            allow_pileup_vcf_fallback=args.allow_pileup_vcf_fallback,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            summary_tsv=args.summary_tsv,
            metadata_json=args.metadata_json,
            extra_clair3_args=args.extra_clair3_arg,
        )
        print(run_clair3_longread_generation(config))
        return


    if args.command == "qc-longread-clair3":
        config = Clair3QcConfig(
            cohort=args.cohort,
            output_root=args.output_root,
            platform=args.platform,
            model_path=args.model_path,
            reference_fasta=args.reference_fasta,
            sample_list=args.sample_list,
            samples=args.sample,
            samtools_bin=args.samtools_bin,
            inspect_reads=args.inspect_reads,
            min_mq=args.min_mq,
            min_coverage=args.min_coverage,
            qual=args.qual,
            min_hp_fraction=args.min_hp_fraction,
            output_prefix=args.output_prefix,
            summary_tsv=args.summary_tsv,
            metadata_json=args.metadata_json,
        )
        result = run_clair3_qc(config)
        print(result)
        # QC JSON/TSV는 먼저 남긴다. 그러나 FAIL 또는 MISSING 상태를 성공으로
        # 반환하면 sample runner가 Stage 06–10으로 진행해 잘못된 VCF/HP 산출물을
        # 사용하게 된다. per-sample runner는 이 non-zero를 해당 sample 실패로만
        # 기록하고 다음 sample으로 계속한다.
        if int(result.get("fail_count", 0)) or int(result.get("missing_count", 0)):
            raise SystemExit(2)
        return



    if args.command == "run-denovo-assembly":
        config = DenovoAssemblyConfig(
            cohort=args.cohort,
            output_root=args.output_root,
            reference_fasta=args.reference_fasta,
            platform=args.platform,
            sample_list=args.sample_list,
            samples=args.sample,
            raw_root=args.raw_root,
            alignment=args.alignment,
            alignment_pattern=args.alignment_pattern,
            contigs=[args.contigs],
            label=args.label,
            assembler=args.assembler,
            hifiasm_bin=args.hifiasm_bin,
            samtools_bin=args.samtools_bin,
            threads=args.threads,
            fastq_threads=args.fastq_threads,
            exclude_flags=args.exclude_flags,
            min_contig_len=args.min_contig_len,
            keep_fastq=args.keep_fastq,
            keep_gfa=args.keep_gfa,
            keep_bin=args.keep_bin,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            summary_tsv=args.summary_tsv,
            metadata_json=args.metadata_json,
        )
        print(run_denovo_assembly(config))
        return

    if args.command == "anchor-denovo-assembly":
        config = AssemblyAnchorConfig(
            cohort=args.cohort,
            output_root=args.output_root,
            reference_fasta=args.reference_fasta,
            sample=args.sample,
            assembly_dir=args.assembly_dir,
            label=args.label,
            minimap2_bin=args.minimap2_bin,
            preset=args.preset,
            threads=args.threads,
            min_mapq=args.min_mapq,
            min_contig_len=args.min_contig_len,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            summary_tsv=args.summary_tsv,
            metadata_json=args.metadata_json,
        )
        print(run_assembly_anchor(config))
        return

    if args.command == "phase-longread-bamonly":
        config = BamOnlyPhaseConfig(
            sample_id=args.sample_id,
            alignment=args.alignment,
            output_dir=args.output_dir,
            reference_fasta=args.reference_fasta,
            contigs=args.contigs,
            platform=args.platform,
            min_mapq=args.min_mapq,
            min_base_quality=args.min_base_quality,
            min_depth=args.min_depth,
            min_marker_reads=args.min_marker_reads,
            min_minor_allele_fraction=args.min_minor_allele_fraction,
            max_major_allele_fraction=args.max_major_allele_fraction,
            max_third_allele_fraction=args.max_third_allele_fraction,
            max_marker_gap=args.max_marker_gap,
            max_markers_per_block=args.max_markers_per_block,
            min_read_markers=args.min_read_markers,
            min_score_delta=args.min_score_delta,
            max_iterations=args.max_iterations,
            include_supplementary=args.include_supplementary,
            exclude_duplicates=args.exclude_duplicates,
            max_pileup_depth=args.max_pileup_depth,
            marker_detection_chunk_bp=args.marker_detection_chunk_bp,
            max_markers_per_contig=args.max_markers_per_contig,
            regions=args.regions,
            unknown_read_policy=args.unknown_read_policy,
            dry_run=args.dry_run,
        )
        print(run_bamonly_phase(config))
        return


    if args.command == "create-personal-realigned-alignments":
        config = PersonalRealignConfig(
            cohort=args.cohort,
            output_root=args.output_root,
            sample_list=args.sample_list,
            samples=args.sample,
            chroms=[args.chroms],
            sequence_variant_policy=args.sequence_variant_policy,
            sequence_source=args.sequence_source,
            assembly_phase_labels_validated=args.assembly_phase_labels_validated,
            dataset_role=args.dataset_role,
            alignment=args.alignment,
            alignment_pattern=args.alignment_pattern,
            hap1_fasta_dir=args.hap1_fasta_dir,
            hap2_fasta_dir=args.hap2_fasta_dir,
            output_prefix=args.output_prefix,
            realigned_subdir=args.realigned_subdir,
            reference_fasta=args.reference_fasta,
            output_format=args.output_format,
            platform=args.platform,
            preset=args.preset,
            threads=args.threads,
            sort_threads=args.sort_threads,
            sort_mem=args.sort_mem,
            jobs=args.jobs,
            hap_jobs=args.hap_jobs,
            route_jobs=args.route_jobs,
            min_mapped_frac=args.min_mapped_frac,
            ratio_qc_policy=args.ratio_qc_policy,
            hap_read_mode=args.hap_read_mode,
            read_routing_policy=args.read_routing_policy,
            unknown_read_policy=args.unknown_read_policy,
            include_unmapped=args.include_unmapped,
            verify_counts=args.verify_counts,
            count_tolerance_abs=args.count_tolerance_abs,
            keep_target=args.keep_target,
            allow_non_pretraining=args.allow_non_pretraining,
            samtools_bin=args.samtools_bin,
            minimap2_bin=args.minimap2_bin,
            bwa_mem2_bin=args.bwa_mem2_bin,
            illumina_aligner=args.illumina_aligner,
            bwa_mem2_index_root=args.bwa_mem2_index_root,
            illumina_duplicate_mode=args.illumina_duplicate_mode,
            inspect_reads=args.inspect_reads,
            dry_run=args.dry_run,
            overwrite=args.overwrite,
            summary_tsv=args.summary_tsv,
            metadata_json=args.metadata_json,
        )
        summary = run_personal_realign_generation(config)
        print(summary)
        return

    if args.command == "create-window-index":
        config = WindowIndexConfig(
            cohort=args.cohort,
            output_root=args.output_root,
            sample_list=args.sample_list,
            samples=args.sample,
            chroms=[args.chroms],
            sequence_variant_policy=args.sequence_variant_policy,
            sequence_source=args.sequence_source,
            window_size=args.window_size,
            offset_unit=1 if args.fixed_window_start else args.offset_unit,
            seed=args.seed,
            drop_n=args.drop_n,
            blacklist_bed=args.blacklist_bed,
            disable_reference_blacklist=args.disable_reference_blacklist,
            max_blacklist_frac=args.max_blacklist_frac,
            jobs=args.jobs,
            overwrite=args.overwrite,
            summary_tsv=args.summary_tsv,
            metadata_json=args.metadata_json,
        )
        summary = run_window_index_generation(config)
        print(summary)
        return

    if args.command == "build-igsr-ont-manifest":
        config = IgsrOntManifestConfig(
            raw_root=args.raw_root,
            output_manifest=args.output_manifest,
            dataset_root=args.dataset_root,
            common_vcf_dir=args.common_vcf_dir,
            reference_fasta=args.reference_fasta,
            contigs=args.contigs,
            dataset_role=args.dataset_role,
            read_type=args.read_type,
            project_id=args.project_id,
            project_label=args.project_label,
            platform=args.platform,
            reference_build=args.reference_build,
            sample_status=args.sample_status,
            disease=args.disease,
            storage_policy=args.storage_policy,
            jobs=args.jobs,
            windows_per_chunk=args.windows_per_chunk,
            compact_mode=args.compact_mode,
            read_coordinate_mode=args.read_coordinate_mode,
            hap1_alignment_pattern=args.hap1_alignment_pattern,
            hap2_alignment_pattern=args.hap2_alignment_pattern,
            haplotagged_alignment_pattern=args.haplotagged_alignment_pattern,
            window_index_pattern=args.window_index_pattern,
            validation_level=args.validation_level,
            compact_compressor=args.compact_compressor,
            compact_compression_level=args.compact_compression_level,
            schedule=args.schedule,
            window_size=args.window_size,
            offset_unit=1 if args.fixed_window_start else args.offset_unit,
            seed=args.seed,
            min_mapq=args.min_mapq,
            max_windows=args.max_windows,
            phase_inspect_max_reads=args.phase_inspect_max_reads,
            whatshap_bin=args.whatshap_bin,
            samtools_bin=args.samtools_bin,
            bcftools_bin=args.bcftools_bin,
            common_vcf_pattern=args.common_vcf_pattern,
            common_vcf_filter_expr=args.common_vcf_filter_expr,
            consensus_filter_expr=args.consensus_filter_expr,
            haplotag_vcf_types=args.haplotag_vcf_types,
            sequence_variant_policy=args.sequence_variant_policy,
            sequence_reference_fasta=args.sequence_reference_fasta,
            synthesize_missing_ps=args.synthesize_missing_ps,
            allow_longread_vcf_haplotag=False,
            allow_unphased_longread=args.allow_unphased_longread,
            samples=args.sample,
            sample_file=args.sample_file,
            pilot_size=args.pilot_size,
            limit=args.limit,
            summary_tsv=args.summary_tsv,
            overwrite_manifest=args.overwrite_manifest,
        )
        summary = build_igsr_ont_manifest(config)
        print(summary)
        return

    if args.command == "preprocess-full-shared-hap":
        config = FullPreprocessConfig(
            sample_id=args.sample_id,
            source_sample_id=args.source_sample_id,
            dataset_role=args.dataset_role,
            read_type=args.read_type,
            project_id=args.project_id,
            project_label=args.project_label,
            platform=args.platform,
            reference_build=args.reference_build,
            sample_status=args.sample_status,
            disease=args.disease,
            hap1_fasta_dir=args.hap1_fasta_dir,
            hap1_chain_dir=args.hap1_chain_dir,
            hap2_fasta_dir=args.hap2_fasta_dir,
            hap2_chain_dir=args.hap2_chain_dir,
            alignment=args.alignment,
            contigs=[item.strip() for item in args.contigs.split(",") if item.strip()],
            dataset_root=args.dataset_root,
            output_sample_root=args.output_sample_root,
            work_root=args.work_root,
            window_index_dir=args.window_index_dir,
            read_coordinate_mode=args.read_coordinate_mode,
            hap1_alignment=args.hap1_alignment,
            hap2_alignment=args.hap2_alignment,
            layout_name=args.layout_name,
            hap_read_mode=args.hap_read_mode,
            storage_policy=args.storage_policy,
            jobs=args.jobs,
            windows_per_chunk=args.windows_per_chunk,
            compact_mode=args.compact_mode,
            validation_level=args.validation_level,
            compact_compressor=args.compact_compressor,
            compact_compression_level=args.compact_compression_level,
            schedule=args.schedule,
            window_size=args.window_size,
            offset_unit=1 if args.fixed_window_start else args.offset_unit,
            seed=args.seed,
            stage70_filter_preset=args.stage70_filter_preset,
            min_mapq=args.min_mapq,
            max_windows=args.max_windows,
            overwrite=args.overwrite,
            resume=not args.no_resume,
            allow_untagged=args.allow_untagged,
            haplotag_if_missing=args.haplotag_if_missing,
            allow_longread_vcf_haplotag=False,
            allow_unphased_longread=args.allow_unphased_longread,
            auto_prepare_hap_fasta=args.auto_prepare_hap_fasta,
            phased_vcf=args.phased_vcf,
            common_vcf_dir=args.common_vcf_dir,
            common_vcf_pattern=args.common_vcf_pattern,
            common_vcf_filter_expr=args.common_vcf_filter_expr,
            consensus_filter_expr=args.consensus_filter_expr,
            haplotag_vcf_types=args.haplotag_vcf_types,
            sequence_variant_policy=args.sequence_variant_policy,
            reference_fasta=args.reference_fasta,
            sequence_reference_fasta=args.sequence_reference_fasta,
            personal_fasta_qc=args.personal_fasta_qc,
            phase_inspect_max_reads=args.phase_inspect_max_reads,
            synthesize_missing_ps=args.synthesize_missing_ps,
            whatshap_bin=args.whatshap_bin,
            samtools_bin=args.samtools_bin,
            bcftools_bin=args.bcftools_bin,
            source_wgs=args.source_wgs,
            source_vcfs=args.source_vcf,
            rna_bigwigs=args.rna_bigwig,
        )
        run_full_preprocess(config)
        return

    if args.command == "validate-compact-v2":
        config = CompactValidationConfig(
            compact_contig_dir=args.compact_contig_dir,
            read_coordinate_mode=args.read_coordinate_mode,
            validation_level=args.validation_level,
            output_json=args.output_json,
        )
        print(run_compact_validation(config))
        return



    if args.command == "validate-stage70-trace":
        contigs = [item.strip() for item in args.contigs.split(",") if item.strip()]
        result = validate_stage70_trace_root(
            args.compact_root, contigs, args.report_root, args.max_examples, args.scratch_dir,
            workers=args.workers,
        )
        print(json.dumps(result["summary"], ensure_ascii=False, sort_keys=True))
        if result["summary"]["status"] != "pass":
            raise SystemExit(2)
        return
    if args.command == "preprocess-full-manifest":
        config = FullManifestConfig(
            manifest=args.manifest,
            status_tsv=args.status_tsv,
            dataset_root=args.dataset_root,
            sample_jobs=args.sample_jobs,
            jobs=args.jobs,
            overwrite=args.overwrite,
            resume=not args.no_resume,
            storage_policy=args.storage_policy,
            sample_catalog=args.sample_catalog,
        )
        run_full_manifest(config)
        return
    parser.error(f"unknown command: {args.command}")


if __name__ == "__main__":
    main()
