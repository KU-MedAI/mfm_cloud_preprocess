# 현재 제한사항

## production을 막는 항목

PacBio file freeze와 raw long-read Stage 04 executor는 완료됐다. 현재 남은 package-level 제한은 다음과 같다.

1. Illumina CRAM/CRAI와 GEUVADIS BAM/BAI는 URL·주 파일 크기는 있으나 package 내 공식 checksum이 없다.
2. HPRC R1 VCF는 graph-derived 표현이므로 donor subset·normalization·REF·phased heterozygous GT QC를 통과하기 전까지 production small-variant source로 확정하지 않는다.
3. HGSVC2 VCF는 공개 TBI가 없어 다운로드 checksum 검증 뒤 `bcftools index -t`로 생성한다.
4. RNA ASE label은 Stage 35 QC의 유전자별 최소 informative count와 적격성 판정을 통과한 항목만 학습 target으로 사용한다.

## 집계 해석

- 기존 RNA URL manifest에는 GEUVADIS BAM URL 462명이 있으나, 그중 13명은 현재 1KGP 3,202 panel 및 최종 3,212명 cohort 밖이므로 strict pair에서 제외했다.
- strict RNA BAM pair 449명은 최종 cohort 안에서 실제 BAM URL이 확인된 집합이다.
- WGS-only는 3,212 - 449 = 2,763명이다.
- phase별 용량은 PacBio read/index·VCF·reference의 동결 byte를 사용한다. Illumina/RNA index size가 없는 부분은 여전히 실제 전송량과 다를 수 있다.

## ASE 제한

RNA BAM의 전체 read가 H1/H2 informative한 것은 아니다. heterozygous marker를 덮지 않는 read는 haplotype을 결정할 수 없다. ASE 모델에는 최소한 informative read 수, H1/H2 count, unassigned 비율, total expression depth, mapping bias와 phase 신뢰도를 함께 제공해야 한다.

## historical package

2026-08-27의 3,202명 패키지는 ONT를 포함한 당시 정책을 설명한다. 현재 package와 sample 수·플랫폼·VCF 정책을 섞어 사용하면 안 된다.

