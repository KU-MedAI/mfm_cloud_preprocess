# 변경 이력

## 2026-09-16: file freeze 및 실행 계약 v2

- PacBio 295명 read 511개, 공개 VCF 관련 파일 704개, exact reference 3개를 URL·bytes·checksum까지 동결
- HG03125를 PacBio CLR로 정정하고 Stage 04/personal realignment에 `map-pb` 분기 추가
- HGSVC2 공개 TBI 부재를 기록하고 다운로드 검증 후 로컬 index 생성
- mixed phase 목록에 PacBio read·project VCF·exact reference를 함께 포함
- MD5, 단일-part S3 ETag, 64 MiB multipart ETag offline 검증 추가
- package schema를 `mfm_share_package_3212_20260916_v2`로 갱신

## 2026-09-16: 3,212명 current package v1

### cohort

- 1KGP 3,202명만을 모수로 한 구성을 PacBio 우선 통합 cohort 3,212명으로 변경
- PacBio 295명 중 285명은 1KGP와 중복, 10명은 1KGP panel 밖
- Illumina-only 2,917명
- ONT production 배정 제거
- HPRC R3 345명 제외 사유와 대체 source 보존

### RNA

- strict GEUVADIS RNA BAM pair 449명 유지
- 새 WGS 선택 기준으로 PacBio 33명, Illumina 416명
- WGS-only 수를 2,763명으로 갱신

### VCF·reference

- 공용 1KGP panel을 모든 long-read에 적용하던 설명 제거
- 프로젝트별 public VCF와 exact reference preflight 도입
- 검증 실패 시 Clair3 자동 fallback 금지

### 준비도

- sample registry와 file manifest를 분리
- PacBio 295명을 `not_frozen`으로 명시
- Illumina/RNA checksum 미동결 상태를 명시
- raw Stage 04 executor의 구현 경계를 명시

### predecessor

- 2026-08-27 `3202_sample_share_package`는 수정하지 않고 역사 자료로 보존

