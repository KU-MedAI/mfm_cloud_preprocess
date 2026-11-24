from __future__ import annotations

import argparse
from pathlib import Path

from .direct_compact_writer import (
    DEFAULT_COMPACT_COMPRESSION_LEVEL,
    DEFAULT_COMPACT_COMPRESSOR,
    SUPPORTED_COMPACT_COMPRESSORS,
)
from .full_preprocess import FULL_LAYOUT_NAME
from .window_index import DEFAULT_BLACKLIST_BED

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
    sample_vcf.add_argument("--extraction-mode", choices=["chromosome", "concat"], default="chromosome", help="개인 변이를 염색체별로 먼저 추출함; concat은 구형 경로")
    sample_vcf.add_argument("--concat-dir", type=Path, help="Directory for <cohort>.autosomes.vcf.gz. Default: --common-vcf-dir")
    sample_vcf.add_argument("--bcftools-bin", default="bcftools")
    sample_vcf.add_argument("--jobs", type=int, default=8, help="Sample-level parallel jobs")
    sample_vcf.add_argument("--concat-threads", type=int, default=4)
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

    assembly_direct = subparsers.add_parser(
        "prepare-provided-assembly",
        help="Validate chromosome-level phased assemblies and adapt them to HAP1/HAP2 chrN FASTA inputs",
    )
    assembly_direct.add_argument("--cohort", required=True)
    assembly_direct.add_argument("--sample", required=True)
    assembly_direct.add_argument("--output-root", required=True, type=Path)
    assembly_direct.add_argument("--reference-fasta", required=True, type=Path, help="GRCh38 provenance/coordinate reference; not used to generate assembly sequence")
    assembly_direct.add_argument("--hap1-fasta", required=True, type=Path)
    assembly_direct.add_argument("--hap1-chrom-alias", required=True, type=Path)
    assembly_direct.add_argument("--hap1-chain", type=Path, help="Official assembly->GRCh38 multi-chain retained as provenance only")
    assembly_direct.add_argument("--hap2-fasta", required=True, type=Path)
    assembly_direct.add_argument("--hap2-chrom-alias", required=True, type=Path)
    assembly_direct.add_argument("--hap2-chain", type=Path, help="Official assembly->GRCh38 multi-chain retained as provenance only")
    assembly_direct.add_argument("--chroms", default="chr1-chr22")
    assembly_direct.add_argument("--phase-label-status", choices=["unverified", "validated"], default="unverified")
    assembly_direct.add_argument("--overwrite", action="store_true")

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
    personal_realign.add_argument("--assembly-phase-labels-validated", action="store_true", help="Allow HP-based routing only after assembly HAP1/HAP2 vs VCF HP1/HP2 orientation validation")
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
    personal_realign.add_argument("--min-mapped-frac", type=float, default=0.90, help="Minimum routed primary mapped fraction")
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
    window_index.add_argument("--disable-reference-blacklist", action="store_true", help="Required for assembly_direct until official multi-chain projection is implemented")
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
    full.add_argument("--sequence-source", choices=["vcf_consensus", "assembly_direct"], default="vcf_consensus")
    full.add_argument("--assembly-phase-labels-validated", action="store_true")
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
