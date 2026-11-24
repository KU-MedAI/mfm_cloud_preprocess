# 다운로드 가이드

## 1. 원칙

sample-level 표의 URL을 임의로 조합하지 않는다. 실제 다운로드는 `manifests/`의 file-level 행만 사용하며 `status`를 반드시 확인한다.

## 2. 현재 가능한 목록

### Illumina 30x 2,918명

`manifests/illumina_30x_files.tsv`에 CRAM과 예측 가능한 CRAI URL이 있다. checksum은 아직 package에 동결되지 않았으므로 원 release manifest에서 MD5를 보완한 뒤 대량 다운로드한다.

```bash
PLATFORM_FILTER=Illumina_30x INCLUDE_RNA=0 \
  bash scripts/01_make_phase_url_lists.sh all

DRY_RUN=1 bash scripts/02_download_phase.sh 1
```

### GEUVADIS RNA BAM 449명

`manifests/rna_geuvadis_files.tsv`에 BAM/BAI URL이 있다. WGS와 같이 받을 때 `INCLUDE_RNA=1`을 사용한다.

```bash
PLATFORM_FILTER=Illumina_30x INCLUDE_RNA=1 \
  bash scripts/01_make_phase_url_lists.sh 1
```

### active PacBio HiFi 294명

`manifests/pacbio_file_manifest.tsv`는 historical 295명/read 511개 freeze를 보존한다. active `pacbio_294_samples.tsv`와 `manifests/pacbio_active_file_manifest.tsv`에는 HiFi 294명/507개 파일만 있으며 downloader가 이 목록과 phase assignment로 filter하므로 HG03125 CLR 4개는 받지 않는다. mixed 또는 PacBio 전용 목록에는 active project의 VCF/TBI/phase tag와 exact FASTA/FAI도 자동 포함된다.

```bash
PLATFORM_FILTER=PacBio INCLUDE_RNA=0 \
  bash scripts/01_make_phase_url_lists.sh all

# 전체 3,212명 mixed phase
PLATFORM_FILTER=all INCLUDE_RNA=1 \
  bash scripts/01_make_phase_url_lists.sh 1
```

## 3. 1KGP Illumina VCF

`scripts/04_download_illumina_vcf.sh`는 Illumina donor용 2022 high-coverage panel만 다룬다. PacBio donor의 프로젝트별 VCF를 대신하지 않는다.

```bash
DRY_RUN=1 CHROMS='1 2' bash scripts/04_download_illumina_vcf.sh
```

## 4. 다운로드 후 검증

1. `bash scripts/03_verify_phase.sh <phase>`로 동결 size와 MD5/S3 multipart ETag 검사
2. CRAM/BAM index 존재 확인과 `samtools quickcheck`
3. 공개 index가 없는 HGSVC2 VCF만 `bcftools index -t`로 TBI 생성
4. header의 sample, read group, contig, M5 확인
5. CRAM decode reference로 실제 decode 시험
6. VCF donor, REF, phase 검사

URL 접근 성공만으로 production 검증을 통과한 것으로 보지 않는다.
