# process_code/v2 — 현재 파이프라인 모식도

```text
manifest + sample list
  └─ run_pipeline.sh / run_sample.sh
      ├─ WGS 01→02→03→[04–05 opt-in]→06→[07–08 opt-in]→09→10→11→12
      │    ├─ WGS-only: 23 validation → 24 cleanup
      │    └─ WGS–RNA pair: WGS cleanup을 미루고 RNA로 진행
      └─ RNA 13→14→15→16→17→18→19→20→21→22→23→24
           → 25 WGS+RNA validation → 26 cleanup

모든 sample 처리 뒤 27_make_catalog.py가 PASS/incomplete/cohort를 집계한다.
```

## WGS 01–12

```text
01 입력 표준화
  BAM/uBAM/CRAM/FASTQ → canonical BAM
       │
02 WGS+RNA preflight
  sample/RG + 22 autosome SQ/LN/M5 + reference/VCF + RNA pair contract
       │
03 VCF 준비 ── [04–05 Clair3 opt-in]
       │
06 reference 좌표 HP/PS haplotag ── [07–08 de novo opt-in]
       │
09 HAP1/HAP2 FASTA·FAI + GRCh38→HAP chain
       │
10 WGS read를 HAP1/HAP2에 재정렬
  Illumina=BWA-MEM2, HiFi/ONT=minimap2
       │
11 HAP별 deterministic 1 MiB windows.tsv
       │
12 Compact Zarr + full validation + 독립 source trace
```

WGS 모델 입력은 `compact/<layout>/<contig>/`의 `hap_sequence.zarr`, 두 `readset` Zarr와 `depth.zarr`다. long-read supplementary evidence는 해당하는 경우 별도 `supplementary_depth.zarr`로 저장한다.

## RNA 13–26

```text
13 RNA BAM/CRAM → paired FASTQ + pair/singleton/discard QC
       │
14 GENCODE GTF → HAP1/HAP2 좌표 GTF·exon + lift QC
       │
15 STAR index → 16 STAR align → 17 Picard mark-only
       │
18 한 번의 depth 계산
  ├─ HAP1/2 coverage.zarr          : 개인 HAP 전장
  └─ HAP1/2 coverage_windows.zarr  : WGS kept=1 window와 같은 순서·좌표
       │
19 RNA realign/QC summary
       │
20 HAP FASTA variant를 WGS read로 검증
       │
21 validated variant를 덮는 fragment의 HAP1/HAP2 상대 AS count
       │
22 HAP별 WGS gene rd/CNV ─────┐
                               ├─ 23 ASE QC
20→21 HAP variant SHA-256 ────┘   continuous/binary·split·mmap·overdispersion
                                      │
24 sample metadata stamp
       │
25 WGS+RNA package validation
       │
26 fail-closed cleanup
```

전장 RNA coverage는 GTEx-compatible 정책으로 duplicate를 포함한다. ASE fragment count와 CNV depth는 QC 목적에 따라 duplicate를 제외한다.

## 검증과 cleanup

```text
WGS-only 23 / RNA pair 25
  ├─ sample ID·mode·metadata
  ├─ WGS window/compact/full/source trace
  ├─ RNA 전장/window Zarr schema·좌표·순서·windows SHA-256
  ├─ HAP variant QC → allelic count meta v5 → ASE QC v3 hash
  └─ 필수 QC·log·provenance
             │
       PASS + failed_checks=[]
             │
WGS-only 24 / RNA pair 26
  삭제 직전에 같은 계약을 현재 파일로 다시 검사
  stale PASS·누락·schema/hash 불일치면 삭제하지 않음
```

최종 보존 대상은 Compact Zarr, RNA HAP별 두 coverage Zarr, HAP variant/ASE 결과, QC, metadata, logs다. 원본 입력은 cleanup 대상이 아니다.

## v2와 TGIL 번호 대응

| 기능 | v2 | TGIL |
|---|---:|---:|
| RNA FASTQ/GTF | 13–14 | 11–12 |
| STAR index/align/MarkDuplicates | 15–17 | 13–15 |
| 전장+WGS-window RNA coverage | 18 | 16 |
| RNA QC | 19 | 17 |
| HAP variants/allelic count/CNV/ASE | 20–23 | 18–21 |
| metadata/validation/cleanup | 24–26 | 22–24 |
| catalog | 27 | 25 |

두 패키지는 번호만 다르고 최근 ASE 계산과 Zarr·검증 계약은 동일하다.
