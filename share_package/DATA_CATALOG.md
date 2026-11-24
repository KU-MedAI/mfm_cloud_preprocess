# 데이터 카탈로그

## 프로젝트별 선택 결과

| 프로젝트 | 선택 인원 | WGS read | 공개 small-variant 자료 | 현재 처리 원칙 |
|---|---:|---|---|---|
| HPRC R2 | 232 | HiFi uBAM, DeepConsensus FASTQ, 공개 aligned BAM | longcallD phased VCF, DeepVariant VCF | exact reference 확인 후 aligned BAM 우선; VCF donor subset·정규화·phase 검증 |
| HGSVC3 | 60 | HiFi uBAM/BAM 또는 archive FASTQ | PAV phased SNV/indel | 모든 run 병합·동일 reference 정렬; GT phase 검증; FORMAT/PS 부재를 별도 기록 |
| HGSVC2 | active 1 | NA12878 HiFi aligned BAM+BAI; HG03125 CLR은 historical freeze만 보존 | PAV Freeze4 phased SNV/indel | 공개 TBI는 checksum 확인 뒤 생성; donor subset·REF/M5 검사·haplotag |
| HPRC R1 | 1 | HiFi CCS 또는 공개 aligned BAM | Minigraph-Cactus graph-derived VCF | HG03492만 선택; graph INFO를 FORMAT/PS로 오해하지 않음 |
| 1KGP Illumina 30x | 2,918 | GRCh38 aligned CRAM+CRAI | 3,202명 phased SNV/indel/SV panel | donor subset에서 SNP/indel 제한; exact CRAM decode reference 검증 |

HPRC R3 345명과 HGSVC1은 현재 production source로 선택하지 않는다.

## read 형식별 Stage 01 진입

| 입력 | 현재 상태 | 필요한 처리 |
|---|---|---|
| aligned CRAM+CRAI | 좌표 정렬 완료 | decode reference와 `@SQ/M5` 확인 후 BAM 정규화 또는 직접 사용 |
| aligned BAM+index | 좌표 정렬 완료 | sort/index, sample/read group, reference 일치 검사 |
| uBAM/CCS BAM | sequence container, `@SQ`가 없을 수 있음 | Stage 01에서 FASTQ stream으로 추출해 HiFi `map-hifi`로 analysis reference에 정렬 |
| FASTQ.gz | raw read | Stage 01에서 Illumina는 BWA-MEM2, HiFi는 minimap2 `map-hifi`로 정렬 후 sort/index. CLR/ONT는 코드가 지원하지만 현재 manifest에는 없음 |

`@SQ=0`은 read sequence가 없다는 뜻이 아니라 reference sequence dictionary가 없어 unaligned container라는 뜻이다.

## VCF 해석

- BAM의 base quality와 VCF의 site `QUAL`은 서로 다른 값이다.
- HGSVC2/3 assembly-derived VCF의 site `QUAL=.`은 variant confidence 숫자가 기록되지 않았다는 뜻이며 read의 base quality 부재를 뜻하지 않는다.
- phased GT가 있으면 haplotype allele 방향은 사용할 수 있지만 FORMAT/PS가 없으면 연속 phase block 경계를 별도 방식으로 다뤄야 한다.
- WhatsHap haplotag는 donor의 phased heterozygous variant와 read가 필요하다. 단순히 VCF 파일이 존재하는 것만으로 충분하지 않다.
- public VCF는 donor 존재, reference REF 일치, phased heterozygous GT, normalization/index를 모두 통과해야 한다.

## RNA

- 기본 RNA pair는 GEUVADIS TOPMed processed GRCh38 BAM이 확인된 449명이다.
- RNA BAM은 expression label, ASE count와 QC에 사용한다.
- recount3 BigWig은 원 read BAM이 아니므로 기본 pair 및 다운로드 대상에서 제외한다.
- MAGE FASTQ, HPRC Kinnex, Iso-Seq은 별도 연구 자료로 취급하며 이 package의 strict RNA pair에 섞지 않는다.
- RNA read 전체가 H1/H2로 배정되는 것은 아니다. informative heterozygous site를 덮는 read만 allele-specific evidence를 제공하므로 ASE는 informative depth, imbalance 검정력, mapping bias를 함께 기록해야 한다.
