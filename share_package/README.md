# MFM WGS/RNA current share package: 3,212명

## 목적

이 디렉터리는 MFM 데이터 후보 목록과 다운로드·전처리 계약을 외부 연구팀에 전달하기 위한 current package다. 2026-09-29 active 선택에서는 CLR을 제외하며, 과거 2026-09-16 PacBio freeze 파일은 provenance 재현용으로 보존한다.

과거 3,202명 패키지와 달리 ONT를 production long-read source로 선택하지 않는다. HPRC R3는 데이터 품질 문제가 아니라 release 변동성과 preprint 단계 재현성 문제로 이번 freeze에서 제외한다.

## 확정 cohort

| 구분 | 인원 |
|---|---:|
| 1KGP Illumina 30x 원 panel | 3,202 |
| 선택 PacBio HiFi | 294 |
| 선택 PacBio CLR | 0 |
| PacBio와 1KGP 중복 | 285 |
| 1KGP 밖 PacBio | 10 |
| 선택 Illumina 30x | 2,918 |
| 최종 고유 개인 | **3,212** |
| strict GEUVADIS RNA BAM pair | **449** |
| WGS-only | **2,763** |

RNA pair 449명의 현재 WGS 선택은 PacBio HiFi 33명, Illumina 30x 416명이다. HG03125는 active cohort에서 동일 개인의 1KGP Illumina 30× CRAM으로 교체했다.

## 가장 중요한 준비도

active PacBio HiFi 294명은 read 파일 507개를 사용한다. `manifests/pacbio_active_file_manifest.tsv`는 이 active HiFi-only 파일 목록이고, 과거 PacBio 295명/read 511개(`pacbio_file_manifest.tsv`) freeze와 공개 VCF 704개·reference 3개는 provenance로 그대로 보존한다. `scripts/00_validate_package.sh`는 active 294/507과 historical 295/511을 각각 검사한다.

- HPRC R2 232명: GRCh38_no_alt aligned BAM+BAI와 개인별 longcallD VCF/TBI/phase TSV.
- HGSVC3 60명: HiFi uBAM/FASTQ와 공유 SNV·indel VCF/TBI.
- HGSVC2 active 1명: NA12878 HiFi aligned BAM+BAI. HG03125 CLR 4개는 historical freeze에만 남고 active 입력은 1KGP Illumina CRAM이다.
- HPRC R1 1명: HG03492 HiFi FASTQ/uBAM과 graph-derived VCF/TBI. 일반 read-based small variant와 동일시하지 않고 normalization/REF/phase QC를 통과해야 한다.
- Illumina 2,918명 CRAM/CRAI와 GEUVADIS RNA 449명 BAM/BAI는 URL이 있지만 package 내 공식 checksum은 아직 미동결이다.

공개 VCF 검증 실패를 Clair3로 조용히 대체하지 않는다.

## 파일 지도

| 파일 | 내용 |
|---|---|
| `samples_3212.tsv`, `.csv` | 최종 3,212명 sample registry와 RNA·준비도 결합표 |
| `pacbio_294_samples.tsv` | active 선택 PacBio HiFi 294명 |
| `illumina_2918_samples.tsv` | active 선택 Illumina 2,918명 |
| `pacbio_295_samples.tsv`, `illumina_2917_samples.tsv` | 2026-09-16 historical selection 기록 |
| `rna_pair_449.tsv` | strict GEUVADIS RNA BAM pair 449명 |
| `hprc_r3_excluded_345.tsv` | HPRC R3 제외 및 대체 source 기록 |
| `project_data_formats.tsv` | 프로젝트별 read/VCF 형식과 전처리 |
| `summary_counts.tsv`, `summary.json` | 주요 집계와 준비도 |
| `phase_assignments.tsv`, `phase_summary.tsv` | 5개 논리 phase와 추정 용량 |
| `manifests/illumina_30x_files.tsv` | 선택된 Illumina CRAM/CRAI URL |
| `manifests/rna_geuvadis_files.tsv` | strict RNA BAM/BAI URL |
| `manifests/pacbio_active_file_manifest.tsv` | active PacBio HiFi 294명 read 507개 file-level 입력 목록 (CLR 0) |
| `manifests/pacbio_file_manifest.tsv` | historical PacBio 295명 read 511개 freeze (HG03125 CLR 4개 포함, provenance 전용) |
| `manifests/public_vcf_sources.tsv` | 프로젝트별 VCF callset/reference 설명 |
| `manifests/public_vcf_file_manifest.tsv` | 공개 VCF/TBI/phase tag 704개 파일 동결표 |
| `manifests/reference_file_manifest.tsv` | exact GRCh38_no_alt FASTA/FAI 및 workflow reference 동결표 |
| `manifests/source_snapshots/` | 공식 index·논문 supplement snapshot과 SHA256 |
| `PACBIO_FREEZE_PROVENANCE.md` | 295명 file freeze 근거와 checksum 규칙 |
| `DATA_CATALOG.md` | 프로젝트별 데이터 형식과 선택 정책 |
| `PIPELINE_GUIDE.md` | 입력 분기, VCF, haplotag, cleanup 흐름 |
| `DOWNLOAD_GUIDE.md` | 안전한 다운로드 절차 |
| `DOWNLOAD_POLICY.md` | release·checksum·phase 운영 원칙 |
| `KNOWN_LIMITATIONS.md` | 현재 미완료 사항과 사용 금지 해석 |
| `CHANGELOG.md` | 3,202명 구판 대비 변경점 |

## 권장 사용 순서

1. `summary.json`과 `KNOWN_LIMITATIONS.md`를 읽는다.
2. `samples_3212.tsv`에서 sample과 준비도 상태를 고른다.
3. URL은 sample 표가 아니라 `manifests/`의 file-level 행을 사용한다.
4. `scripts/01_make_phase_url_lists.sh`로 read·VCF·reference 목록을 만들고 다운로드 후 `03_verify_phase.sh`로 동결값을 검증한다.
5. `PIPELINE_GUIDE.md`의 입력 형식 및 public VCF preflight를 통과시킨다.
6. 대표 sample smoke test 후 전체 phase를 실행한다.

## 버전

- package schema: `mfm_share_package_3212_20260929_v3`
- active cohort: 2026-09-29
- historical PacBio resource freeze: 2026-09-16
- historical predecessor: `260826/3202_sample_share_package`

## WGS-RNA pair 우선 전처리

새 pair-first 모드는 strict RNA BAM pair를 먼저 완결한 뒤 WGS-only를 처리한다. 선택 단계에서는 local 후보 TSV의 WGS alignment/인덱스와 RNA BAM/인덱스를 모두 다시 확인하고, pair donor를 먼저 고정해 WGS-only에서 제외한다. RNA 또는 index가 빠진 후보는 WGS-only로 자동 대체하지 않으며 selection QC의 부족분·사유에 남긴다.

phase 목록만 만들려면 MODE=pair_first로 scripts/05_preprocess_phase.sh를 실행한다. 각 phase에는 pair WGS, pair RNA, WGS-only 순서를 기록한 execution_order TSV가 생성된다. 실제 실행은 RUN_PREPROCESS=1과 PAIR_WGS_CMD, PAIR_RNA_CMD, WGS_ONLY_CMD를 모두 명시할 때만 수행된다. 전체 phase에서는 모든 pair WGS, 모든 pair RNA, 모든 WGS-only 순으로 실행한다. 세 명령은 기존 sample-isolated wgs-list/rna-list 경로를 사용해야 하며, sample validation 실패는 해당 runner의 QC에 남기고 다음 sample으로 진행해야 한다.
