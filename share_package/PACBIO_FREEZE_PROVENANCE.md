# PacBio 295명 파일 동결 근거

동결일은 2026-09-16이다. 이 문서는 `pacbio_295_samples.tsv`의 295명을 실제 다운로드 가능한 read·variant·reference 파일에 연결한 기준을 기록한다. HPRC R3는 preprint/release tracker 단계 자료를 첫 production cohort에서 제외한다는 결정에 따라 포함하지 않았다.

## 동결 결과

| 선택 프로젝트 | 개인 수 | read 파일 행 | 실제 read 기술 | 주 입력 |
|---|---:|---:|---|---|
| HPRC R2 | 232 | 232 | PacBio HiFi | GRCh38_no_alt aligned BAM+BAI |
| HGSVC3 | 60 | 267 | PacBio HiFi | uBAM 또는 FASTQ.gz |
| HGSVC2 | 2 | 5 | HG03125 CLR 4개, NA12878 HiFi 1개 | CLR FASTQ.gz 또는 aligned BAM+BAI |
| HPRC R1 | 1 | 7 | PacBio HiFi | DeepConsensus FASTQ.gz 3개+Revio uBAM 4개 |
| 합계 | 295 | 511 | HiFi 507개+CLR 4개 | 모든 행 URL·bytes·checksum ready |

HGSVC3 run은 Nature 2025 논문의 Supplementary Table 2에서 선택 60명의 `use data=yes` HiFi accession만 추출했다. 보충표에 같은 accession이 반복된 경우에는 한 번만 유지했다. 각 accession은 ENA file report의 submitted primary file을 우선 사용하고, submitted file이 없으면 archive FASTQ를 사용했다. `SRR19416256`은 ENA file URL이 비어 있어 보충표의 파일명을 공식 IGSR HG03520 sample manifest와 수동 매칭했다.

HG03125는 HGSVC2 35-genome cohort의 PacBio 자료이지만 HiFi가 아니라 CLR/subread이다. 따라서 전체 플랫폼을 PacBio라고 묶을 수는 있어도 295명 전원을 HiFi라고 부르면 안 된다. Stage 04는 이 네 파일에 `minimap2 -x map-pb`를 사용해야 한다. NA12878은 GIAB의 Q20 CCS 약 30× GRCh38 haplotagged BAM을 사용한다.

## Variant와 reference

- HPRC R2: 232명 각각의 GRCh38_no_alt longcallD VCF, TBI, HP/PS 복원 TSV를 동결했다. 총 696개 파일이다.
- HGSVC3: GRCh38 SNV VCF/TBI와 indel VCF/TBI 4개를 직접 내려받아 MD5를 계산했다.
- HGSVC2: Freeze4 SNV·indel VCF 2개와 공식 manifest MD5를 기록했다. 공개 TBI가 없으므로 다운로드 후 `bcftools index -t`로 생성한다.
- HPRC R1: minigraph-cactus decomposed graph VCF와 TBI를 동결했다. HG03492가 header sample에 존재하지만, 일반 read-based small-variant VCF와 동일하다고 간주해서는 안 된다. donor subset·normalization·REF 검사 후 사용한다.
- 공통 exact reference는 NCBI `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`와 FAI를 MD5까지 고정했다. HPRC R2 workflow의 공식 `GRCh38_no_alt.tar.gz`도 별도로 고정했다.

VCF가 공개되어 있다는 이유만으로 무조건 Clair3를 생략하지 않는다. HPRC R2 longcallD와 HGSVC2/3 SNV·indel callset은 공개 VCF 분기를 사용할 수 있다. HPRC R1 graph VCF는 sample 추출과 표현 정규화 QC가 통과할 때만 사용하며 실패하면 별도 small-variant calling이 필요하다.

## checksum 규칙

- `md5`: 제공기관 MD5 또는 이 동결 작업에서 전체 파일을 내려받아 계산한 MD5이다.
- `md5_s3_etag`: 단일-part S3 object ETag이므로 MD5로 검증한다.
- `s3_multipart_etag_64MiB`: 64 MiB part 크기로 재계산 가능한 multipart ETag이다. 일반 `md5sum`과 값이 같지 않다.
- HPRC R2 HG002 aligned BAM은 공식 index에 `.bam.md5` 경로가 있으나 object가 404이므로 versioned S3 multipart ETag를 기록했다. 나머지 HPRC R2 BAM 231개는 공식 `.bam.md5`를 사용한다.

## 검증과 갱신

`scripts/06_validate_pacbio_freeze.py`는 네트워크 없이 295명 donor set, 511개 read 파일, 704개 VCF 관련 파일, 3개 reference resource, URL·bytes·checksum 형식, 프로젝트별 기술 구분을 fail-closed로 검사한다. 원본 공식 index와 보충표에서 다시 생성할 때는 archive 계획 폴더의 `build_pacbio_freeze.py`와 `build_public_resources.py`를 사용하고, 기존 manifest를 먼저 백업한다.
