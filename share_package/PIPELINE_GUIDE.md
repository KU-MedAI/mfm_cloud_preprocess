# 현재 전처리 파이프라인 가이드

## 1. 전체 흐름

```text
sample registry
  -> file-level manifest와 checksum freeze
  -> 입력 preflight(reference/read group/index)
  -> Stage 04용 analysis BAM 정규화
       aligned CRAM/BAM -> 검증·정규화
       uBAM/FASTQ       -> analysis reference 정렬
  -> 개인별 public VCF preflight
       통과             -> donor SNP/indel subset·정규화
       VCF 부재         -> 승인된 platform model로 Clair3 검토
       검증 실패        -> 중단(자동 fallback 금지)
  -> WhatsHap haplotag(HP/PS)
  -> H1/H2 personal FASTA·window
  -> WGS personal realignment
  -> compact model input 생성·검증
  -> WGS-RNA pair는 RNA 정렬/ASE/QC 완료
  -> 조건 충족 후 중간 BAM cleanup
```

## 2. 입력 계약

각 sample은 최소한 다음 정보를 고정해야 한다.

- sample/donor ID와 alias
- read 형식, 정렬 여부, 플랫폼
- analysis reference FASTA
- CRAM이면 decode reference
- public VCF 파일/TBI, donor ID, callset와 exact VCF reference
- raw 입력이면 최종 prepared alignment 경로
- WGS-RNA pair이면 RNA BAM/BAI와 RNA reference

현재 v2의 `01_prepare_analysis_bam.sh`는 aligned BAM/CRAM을 QC 후 canonical BAM으로 표준화한다. raw Illumina paired FASTQ/uBAM은 BWA-MEM2 `mem -Y`, HiFi/CLR/ONT는 minimap2 `map-hifi/map-pb/map-ont`로 정렬·병합·sort·index한다. `@SQ`, sample `SM`, coordinate sort, index, base quality, mapped read와 record/tag 수를 fail-closed로 검사한다.

Illumina duplicate는 canonical 좌표에서 한 번만 확정한다. 기본 `auto`는 source의 비삭제 Picard MarkDuplicates 이력을 검증해 기존 0x400을 보존하고, raw 또는 이력 없는 입력만 Picard `REMOVE_DUPLICATES=false`로 표시한다. Stage 02는 metrics 또는 source provenance를 확인하며 불일치 sample만 중단한다. HAP별 개인 좌표에서는 재표시하지 않고 source flag를 운반하며 최종 short-read model input에서 제외한다. 표준 HiFi/ONT WGS에는 blanket Picard marking을 적용하지 않는다.

## 3. public VCF와 Clair3

현재 선택 프로젝트는 모두 사용할 수 있는 공개 small-variant 자료원이 있으므로 우선 public VCF 분기를 사용한다. 그러나 project root가 있다는 사실과 특정 donor가 production-ready라는 것은 다르다.

public VCF 분기 통과 조건:

1. donor가 header에 존재
2. BAM/CRAM과 VCF의 reference build 및 contig/M5 일치
3. REF allele 검사 통과
4. phased heterozygous GT 존재
5. SNP/indel 제한, multiallelic 분해와 normalization 완료
6. bgzip/index 정상

실패 시 Clair3로 자동 fallback하지 않는다. VCF가 실제로 없는 경우에만 플랫폼·모델·reference를 명시해 별도 calling을 승인한다.

## 4. haplotag와 RNA/ASE

- WGS read는 phased VCF와 WhatsHap을 이용해 HP/PS를 부여한다.
- RNA read의 haplotype 증거는 informative variant를 덮는 read에 한정된다.
- ASE 결과에는 H1/H2 count뿐 아니라 total depth, informative read 비율, unassigned read, mapping bias QC를 저장한다.
- WGS와 RNA의 sample ID, reference, phase orientation이 일치하지 않으면 pair 처리를 중단한다.

## 5. cleanup 시점

- WGS-only: compact output 검증과 checksum 기록 후 재생성 가능한 중간 BAM을 제거할 수 있다.
- WGS-RNA pair: RNA personal alignment, ASE, post-training compact 검증이 끝날 때까지 WGS haplotagged/personal BAM을 보존한다.
- raw 원본, 공개 VCF, exact reference와 최종 compact output은 자동 cleanup 대상이 아니다.

## 6. 현재 v2 문서

실행 wrapper와 CLI 세부 명령은 서버의 `process_code/v2/README.md`를 기준으로 한다. 이 공유 패키지는 cohort와 입력 계약을 전달하며 서버 절대경로나 실행 중인 job 상태를 포함하지 않는다.

