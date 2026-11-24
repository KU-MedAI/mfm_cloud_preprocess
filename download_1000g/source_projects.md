# Source projects and download links

> **상태:** 아래 표는 historical 3,202명 패키지가 참조한 자료원이다. 현재 선택 인원과 project별 public VCF/reference 준비도는 `../share_package/DATA_CATALOG.md` 및 `../share_package/manifests/public_vcf_sources.tsv`를 사용한다.

이 문서는 3,202 sample 공유 패키지에서 사용하는 원천 project와 다운로드 위치를 정리한 파일임. 모든 package 내부 파일 경로는 상대경로로 작성함.

## Project table

| label | data type | file format | official project/root link | package manifest | note |
|---|---|---|---|---|---|
| 1000 Genomes Project high-coverage 30x WGS / ENA PRJEB31736, PRJEB36890, PRJEB55077 | Illumina short-read WGS | `*.final.cram`, `*.final.cram.crai` | https://www.internationalgenome.org/data-portal/data-collections/1000genomes_30x/ | `manifests/illumina_30x_ena_final_cram_size_by_sample.tsv` | NYGC 30x 3,202 sample 공식 collection임. PRJEB31736은 2,504 panel, PRJEB36890은 698 related sample, PRJEB55077은 3,202/602 trios 통합 accession으로 함께 표기하는 것이 정확함. |
| 1KG_ONT_VIENNA | ONT WGS aligned to hg38 | `*.hg38.cram`, `*.hg38.cram.crai` | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/ | `manifests/wgs_urls.tsv` | 현재 전처리에 바로 넣기 좋은 ONT aligned CRAM source임. |
| Human Pangenome Reference Consortium ONT | ONT raw/read BAM | BAM/uBAM multi-run | https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/ | `manifests/wgs_urls.tsv` | aligned CRAM이 아니라 raw uBAM 성격이면 hg38/GRCh38 alignment가 필요함. |
| Human Pangenome Reference Consortium PacBio HiFi | PacBio HiFi WGS | BAM/uBAM multi-run | https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/ | `manifests/wgs_urls.tsv` | sample별 여러 run이 있을 수 있음. manifest에서 sample/source로 필터링해서 전체 run을 받아야 함. |
| HGSVC3 PacBio HiFi/CCS | PacBio HiFi/CCS WGS | BAM/uBAM multi-run | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/ | `manifests/wgs_urls.tsv` | WGS HiFi/CCS run만 사용함. `Kinnex` 또는 `*.flnc.bam`은 RNA 성격이므로 WGS 다운로드에서 제외함. |
| GEUVADIS TOPMed RNA-seq | processed RNA BAM | STAR GRCh38 BAM | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/geuvadis/working/geuvadis_topmed/data/ | `manifests/rna_urls.tsv` | strict RNA-pair의 주 source임. |
| 1000 Genomes Project 2022 phased panel | phased SNV/INDEL/SV VCF | `*.vcf.gz`, `*.vcf.gz.tbi` | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/ | `README.md`, `samples_3202.csv` | H1/H2 sequence 생성 및 haplotagging marker로 사용함. |

## Download examples using package-relative manifests

### Illumina 30x WGS CRAM

```bash
awk -F'\t' '$1=="HG00096"{print "https://"$3}' manifests/illumina_30x_ena_final_cram_size_by_sample.tsv
```

### ONT Vienna aligned CRAM

```bash
awk -F'\t' '$1=="HG00096" && $2=="ONT_Vienna"{print $6}' manifests/wgs_urls.tsv
```

### HPRC HiFi all runs for one sample

```bash
awk -F'\t' '$1=="HG00099" && $2=="HPRC_HiFi"{print $6}' manifests/wgs_urls.tsv
```

### GEUVADIS RNA BAM

```bash
awk -F'\t' '$1=="HG00096" && $2=="Geuvadis_TOPMed"{print $6}' manifests/rna_urls.tsv
```

## Interpretation notes

- `samples_3202.csv`의 `wgs_primary_url`은 대표 URL임. HiFi/CCS와 HPRC ONT는 multi-run일 수 있으므로 반드시 `wgs_url_count_for_assigned_platform`와 manifest를 확인함.
- strict RNA-pair는 RNA BAM 존재 기준임. BigWig은 이번 다운로드/전처리 대상에서 제외함.
- `missing_or_lookup_required`는 donor matrix에는 platform flag가 있지만 현재 package manifest에는 URL이 없는 경우임.
- package에 포함된 manifest는 연구 당시 local audit 결과임. 실제 대량 다운로드 전에는 HTTP HEAD 또는 checksum으로 remote availability를 다시 검증하는 것이 안전함.

## Officialness notes

- `ftp.1000genomes.ebi.ac.uk`와 `internationalgenome.org`는 IGSR/1000 Genomes 공식 배포 경로임.
- ENA accession은 원시/정렬 sequencing data의 공식 archival accession으로 봄.
- HPRC `human-pangenomics` S3 bucket은 HPRC 공식 공개 cloud bucket으로 봄.
- HGSVC3 `data_collections/HGSVC3`는 IGSR FTP에 올라온 HGSVC 공식/영구 배포 경로이나, `working/` 하위 파일은 publication-ready release와 preliminary/working file이 섞일 수 있어 data reuse statement와 파일명 확인이 필요함.
- recount3 BigWig은 원본 GEUVADIS BAM이 아니라 recount3가 균일하게 재처리한 RNA-seq signal track임. 공식 원천 RNA BAM이 아니며 이번 공유 패키지의 다운로드 대상에서 제외함.
- MAGE RNA FASTQ와 HPRC Kinnex/FLNC RNA는 RNA BAM 기준 운영 대상에서 제외함.

## Excluded resources

| label | reason | previous manifest |
|---|---|---|
| recount3 GEUVADIS/SRA base_sums BigWig | 원본 GEUVADIS RNA BAM이 아니라 재처리 signal track이며 이번 운영에서 다운로드하지 않음. | archive package에만 보존 |
| MAGE RNA FASTQ / HPRC Kinnex RNA | RNA BAM 기준이 아니므로 이번 운영에서 다운로드하지 않음. | archive backup 또는 과거 manifest에만 보존 |
