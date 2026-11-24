# 다운로드·release 동결 정책

## cohort freeze

- package version: `mfm_share_package_3212_20260916_v2`
- PacBio 우선순위: HPRC R2 → HGSVC3 → HPRC R1 → HGSVC2
- ONT와 HPRC R3는 이번 production 선택에서 제외
- 동일 개인은 선택 WGS source 하나만 사용

## file freeze 필수 열

각 실제 파일은 다음 값이 있어야 `ready`가 된다.

- canonical sample ID
- project/release/accession
- 직접 file URL
- file type과 alignment state
- size와 MD5 또는 공식 checksum
- reference build 및 가능한 경우 FASTA checksum
- BAM/CRAM index 또는 VCF TBI

active PacBio HiFi 294명의 모든 선택 파일은 `ready`로 동결됐다. historical manifest에는 HG03125 CLR을 포함한 295명/511개를 보존하지만 active downloader는 해당 4개 CLR 파일을 선택하지 않는다.

## VCF 정책

- Illumina cohort VCF는 Illumina 선택 donor에만 기본 적용한다.
- PacBio는 선택 프로젝트의 donor-specific 또는 multisample callset을 사용한다.
- project root만 확인된 상태를 file-ready로 표시하지 않는다.
- public VCF 검증 실패 시 Clair3로 자동 fallback하지 않는다.
- VCF site `QUAL=.`만으로 read base quality가 없다고 판단하지 않는다.

## phase/shard

- 기본 논리 phase는 5개다.
- Illumina/RNA는 manifest size를 사용한다.
- PacBio는 read/index와 project VCF/reference의 동결된 실제 byte를 사용한다. shared resource는 phase 안에서 URL 중복 제거 후 한 번만 센다.
- Illumina/RNA index처럼 공식 size가 없는 일부 파일 때문에 phase 총량은 여전히 근사치다.

## RNA

- strict pair 449명의 GEUVADIS processed BAM만 기본 포함한다.
- BigWig, MAGE, Kinnex, Iso-Seq은 기본 다운로드에서 제외한다.
- RNA가 있는 sample의 WGS BAM은 RNA/ASE 및 post-training QC가 끝날 때까지 cleanup하지 않는다.
