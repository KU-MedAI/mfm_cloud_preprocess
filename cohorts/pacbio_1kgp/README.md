# PacBio 우선·1KGP WGS 통합 cohort

## 1. 이 디렉터리의 역할

이 디렉터리는 HPRC·HGSVC PacBio와 1KGP Illumina 30× WGS를 개인 단위로 합친 **cohort registry**다. 동일 개인은 한 행만 유지하고 선택 PacBio와 Illumina가 모두 있으면 PacBio를 선택한다.

이 디렉터리의 TSV는 cohort registry다. 2026-09-16에 read·VCF·reference URL·bytes·checksum을 동결한 실행 자료는 `../../share_package/manifests/`에 있으며, 로컬 다운로드 경로를 채운 실행용 `PIPELINE_MANIFEST`는 phase 다운로드 후 생성한다.

2026-09-16 외부 전달용 current snapshot은 `../../share_package/`에 있다. 이 registry에 strict GEUVADIS RNA pair, 다운로드 준비도, project VCF 상태와 fail-closed 스크립트를 결합한 패키지다.

## 2. 최종 cohort

| 구분 | 인원 |
|---|---:|
| 1KGP Illumina 30× 원 panel | 3,202 |
| production 고유 PacBio | 294 |
| └ HiFi / CLR | 294 / 0 |
| PacBio와 1KGP panel 중복 | 284 |
| 1KGP panel 밖 PacBio | 10 |
| PacBio로 대체된 1KGP 개인 | 284 |
| Illumina 선택 | 2,918 |
| **최종 고유 개인** | **3,212** |

계산은 `3,202 - 284 + 294 = 3,212`다. 최종 목록에서 `sample_id` 중복은 0건이다. HG03125는 historical CLR freeze에는 남지만 production 선택은 1KGP Illumina다.

1KGP panel 밖 PacBio 10명도 최종 목록에 포함하며 `in_1kgp_3202_panel=N`으로 구별한다.

`HG002, HG005, HG01123, HG02109, HG02486, HG02559, HG03471, HG06807, NA21309, NA21487`

## 3. 최종 자료원 선택

| 선택 프로젝트 | 인원 | read |
|---|---:|---|
| HPRC R2 | 232 | PacBio HiFi |
| HGSVC3 | 60 | PacBio HiFi |
| HPRC R1 | 1 | PacBio HiFi |
| HGSVC2 | 1 | NA12878 HiFi; HG03125 CLR은 historical freeze만 보존 |
| HGSVC1 | 0 | CLR 중심, 선택 없음 |
| 1KGP Illumina 30× | 2,918 | PCR-free short read |

PacBio 내부 중복 선택 우선순위는 다음과 같다.

`HPRC R2 → HGSVC3 → HPRC R1 → HGSVC2 → HGSVC1`

이 우선순위는 프로젝트 이름의 최신성만이 아니라 HiFi read, aligned BAM, small-variant VCF, phase 정보와 exact reference의 연결 가능성을 기준으로 한다.

## 4. 파일 구성

### `pacbio_1kgp_3212_samples.tsv`

최종 3,212명 전체 registry다. 주요 열은 다음과 같다.

- `sample_id`, `aliases`: canonical donor ID와 alias
- `selected_platform`, `selected_project`, `selected_release`: 최종 선택 자료원
- `in_1kgp_3202_panel`: 실제 1KGP 3,202 panel 포함 여부
- `all_published_pacbio_memberships`: 개인이 속한 안정 PacBio cohort
- `hprc_r3_member_but_not_selected`, `hprc_r3_exclusion_reason`: R3 membership과 제외 근거
- `read_type`, `read_format`, `read_alignment_state`, `read_reference`: read 입력 상태
- `read_primary_url_or_rule`, `read_source_index`: 공식 자료원 또는 파일 선택 규칙
- `small_variant_vcf_available`, `vcf_scope`, `vcf_format`, `vcf_callset`: VCF 제공 상태
- `vcf_phase_state`, `vcf_reference`, `vcf_source`: phase와 reference provenance
- `required_preprocessing`, `pipeline_readiness`: 실행 전 필요한 처리

### `pacbio_294_samples.tsv`와 historical `pacbio_295_samples.tsv`

현재 production PacBio 294명은 모두 HiFi이며 `../../share_package/pacbio_294_samples.tsv`와 `../../share_package/manifests/pacbio_active_file_manifest.tsv`를 사용한다. `pacbio_295_samples.tsv`와 `pacbio_file_manifest.tsv`는 HG03125 CLR을 포함한 historical freeze/provenance로만 보존한다.

### `hprc_r3_excluded_345.tsv`

HPRC R3 345명의 제외 사유와 대신 선택된 안정 자료원을 기록한다. R3 개인 345명은 다음과 같이 모두 대체된다.

- HPRC R2 PacBio: 200명
- 1KGP Illumina: 145명

### `project_data_formats.tsv`

프로젝트별 cohort 수, 최종 선택 수, read 형식, small-variant VCF와 필요한 전처리를 요약한다.

## 5. HPRC R3 제외 정책

HPRC R3는 데이터 품질이 나빠서 제외한 것이 아니다. 공식 저장소가 release tracker 형태이고 cohort manifest와 파일이 갱신될 수 있어 첫 production freeze의 재현성을 보장하기 어렵기 때문이다.

- R3 파일은 다운로드 source로 선택하지 않는다.
- R3 membership은 provenance로 보존한다.
- 이후 고정 release tag, accession, checksum과 논문 상태가 확정되면 별도 version으로 다시 평가한다.
- 기존 3,212명 목록을 덮어쓰지 말고 새 cohort version을 생성한다.

주의: HPRC R2 공식 release 안내도 `pre-publication form`이라는 표현을 사용한다. 다만 sequencing index v1.0이 동결되어 있고 현재 입력 준비도가 가장 높으므로 이번 cohort에는 포함했다. 문서·논문에서는 “출판 완료 cohort”가 아니라 “version이 고정된 공식 pre-publication data release”로 표현하고 index version과 checksum을 고정해야 한다.

## 6. 프로젝트별 파이프라인 입력

### HPRC R2

- read: unaligned HiFi BAM, DeepConsensus FASTQ.gz, GRCh38_no_alt/CHM13v2 aligned BAM
- VCF: 개인별 longcallD phased VCF와 DeepVariant small-variant VCF
- 권장: GRCh38_no_alt aligned BAM과 longcallD VCF를 사용한다.
- longcallD phase TSV로 HP/PS를 복원하거나 검증된 phased VCF로 WhatsHap haplotag를 수행한다.
- unphased DeepVariant VCF만으로 public phased 분기에 넣지 않는다.

### HGSVC3

- read: sample별 다수 PacBio HiFi uBAM/BAM 또는 archive FASTQ.gz
- VCF: 65명 PAV v2.4.0.1 assembly-derived phased SNV/indel multisample VCF
- GT는 phased지만 FORMAT/PS가 없고 site QUAL은 `.`이다.
- sample의 모든 run을 합친 후 exact GRCh38-NoALT reference에 정렬하고 donor VCF를 추출·정규화한다.
- 현재 기존 로컬 URL 목록은 65명 중 34명만 포함하므로 production 다운로드에 사용하지 않는다.

### HPRC R1

- read: HiFi CCS BAM 및 GRCh38/CHM13 Winnowmap aligned BAM
- VCF: HPRC v1.0 Minigraph-Cactus graph-derived VCF
- 최종 선택 donor는 HG03492 한 명이다.
- graph VCF의 INFO/PS를 FORMAT phase-set으로 해석하지 않는다.

### HGSVC2

- 전체 cohort는 PacBio CLR/HiFi 혼합이다.
- 최종 선택 HG03125는 CLR FASTQ 4개를 `map-pb`로 정렬하며, NA12878은 GIAB HiFi GRCh38 aligned BAM+BAI를 사용한다.
- PAV Freeze4 SNV/indel multisample VCF는 phased GT를 가지지만 FORMAT/PS가 없고 site QUAL은 `.`이다.
- 공식 `hg38.no_alt`와 read alignment/VCF REF를 대조한 뒤 donor subset과 WhatsHap haplotag를 수행한다.

### HGSVC1

- PacBio CLR 중심이며 SV/assembly 자료에 초점이 있다.
- 이번 목록에서는 선택 donor가 0명이므로 별도 CLR 재정렬·small-variant calling을 수행하지 않는다.

### 1KGP Illumina 30×

- read: GRCh38 aligned CRAM+CRAI
- VCF: 3,202명 phased SNV/indel/SV cohort panel
- active PacBio와 중복되지 않는 2,918명을 Illumina로 선택한다.
- exact CRAM decode reference를 고정하고 cohort VCF에서 개인별 SNP/indel을 추출한다.

## 7. v2 실행 manifest로 변환하는 방법

v2의 현재 실행 manifest는 `mfm_preprocess_manifest_v3` schema다. 이 registry에서 다음과 같이 변환한다.

| cohort registry | 실행 manifest |
|---|---|
| `sample_id` | `sample_id` |
| PacBio HiFi `selected_platform` | `platform=hifi`, `read_type=long_read` |
| Illumina `selected_platform` | `platform=illumina`, `read_type=short_read` |
| `selected_project` | `project_id`, `cohort` |
| 검증·다운로드한 BAM/CRAM | `alignment_path` |
| 공개 phased VCF | `vcf_source=public`, `public_vcf_path` |
| donor ID | `public_vcf_sample_id` |
| exact VCF reference | `vcf_reference_fasta` |
| `vcf_callset` | `vcf_callset` |
| 제공된 phase | `phase_source=provided` |

변환 전에 반드시 다음 필드를 실제 파일 경로로 확정해야 한다.

1. frozen file manifest에서 phase에 해당하는 read·VCF·reference 다운로드
2. 다운로드·정렬된 `alignment_path` 또는 raw `read_files_manifest`
3. VCF 및 TBI의 로컬 경로
4. 동결된 exact reference FASTA의 로컬 경로
5. output `final_root`, `work_root`

public VCF 분기는 다음 조건을 모두 통과해야 한다.

- donor가 VCF header에 존재
- BAM/CRAM과 VCF reference가 동일
- REF allele 전수 검사 통과
- phased heterozygous GT 존재
- SNP/indel 정규화 및 index 정상

검증 실패 시 public VCF를 조용히 Clair3 결과로 대체하지 않는다. 원인을 기록하고 해당 sample을 중단해야 한다.

## 8. 다음 산출물

file freeze는 `../../share_package/`에 완료됐다. 다음 운영 산출물은 다운로드 경로를 반영한 `pipeline_manifest.preflight.tsv`와 대표 sample preflight 통과 후의 `pipeline_manifest.production.tsv`다.

sample registry와 file manifest를 분리해야 한 sample의 여러 sequencing run을 누락하지 않고, cohort 구성 변경 없이 다운로드 파일만 갱신할 수 있다.

## 9. 공식 자료원

- 1KGP 30×: https://internationalgenome.org/data-portal/data-collections/1000genomes_30x/
- HPRC R1: https://github.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0
- HPRC R2: https://humanpangenome.org/hprc-data-release-2/
- HPRC R2 sequencing index: https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables/sequencing_data
- HPRC R2 variant calls: https://github.com/wwliao/hprc_release2_variant_calling
- HPRC R3 tracker: https://github.com/human-pangenomics/R3_HPRC
- HGSVC2: https://internationalgenome.org/data-portal/data-collections/hgsvc2/
- HGSVC3: https://internationalgenome.org/data-portal/data-collections/hgsvc3/

## 10. 생성 기록

- 생성일: 2026-09-14
- 생성 환경: `conda jj_mfm`
- 원본 기록: `/NHNHOME/korea/pre-training/jungjin/codex_archive/plan/260914/`
- 현재 production 검증값: 최종 3,212명, 중복 0건, PacBio HiFi 294명, 1KGP panel 중복 284명, Illumina 2,918명
- historical freeze: PacBio 295명(HiFi 294, CLR 1), 중복 285명 기록을 provenance로 보존

## 11. Library PCR provenance

실행 manifest는 선택 플랫폼이 아니라 프로젝트 Methods에 근거해 증폭 이력을 기록한다. 1KGP 2,918명은 `pcr_free_confirmed`, HPRC R1/HGSVC2/HGSVC3 62명은 `pcr_free_protocol`, HPRC R2 232명은 sample별 PCR field 부재 때문에 `pcr_free_inferred`다. RNA pair 449명은 GEUVADIS TruSeq RNA v2의 10-cycle PCR을 근거로 `pcr_amplified_confirmed`다.

이 분류는 duplicate filtering과 분리한다. 1KGP upstream Picard flag는 보존하고 PacBio에는 blanket MarkDuplicates를 적용하지 않는다. 전체 근거와 manifest 필드는 [../../LIBRARY_PROVENANCE.md](../../LIBRARY_PROVENANCE.md)를 참조한다.
