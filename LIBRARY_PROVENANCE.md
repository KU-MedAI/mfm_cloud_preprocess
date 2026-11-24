# Library 증폭 이력과 duplicate 처리 정책

## 핵심 원칙

PCR 여부는 sequencing platform이나 duplicate flag로 추정하지 않고 프로젝트의 Methods·공식 protocol로 판정한다. Library 증폭 이력과 duplicate 처리는 서로 다른 정보다. PCR-free library에도 optical duplicate, patterned-flow-cell duplicate, coordinate duplicate가 생길 수 있으며, duplicate flag만으로 wet-lab PCR 여부를 역추론할 수 없다.

## 현재 cohort 판정

| 프로젝트 | WGS 표본 | manifest 상태 | 근거와 해석 |
|---|---:|---|---|
| 1KGP_ILLUMINA | 2,918 | `pcr_free_confirmed` | TruSeq DNA PCR-Free를 Methods가 명시한다. 원 workflow의 Picard MarkDuplicates flag는 계속 보존한다. |
| HPRC_R1 | 1 | `pcr_free_protocol` | HMW DNA, SMRTbell Express 2.0, size selection에 PCR 단계가 없다. |
| HPRC_R2 | 232 | `pcr_free_inferred` | 표준 HiFi WGS는 PCR-free이나 공개 index에 sample별 PCR/kit 열이 없어 확정으로 승격하지 않는다. DeepConsensus는 계산적 consensus이며 PCR이 아니다. |
| HGSVC2 | 1 | `pcr_free_protocol` | 선택된 NA12878의 SMRTbell/HiFi Express와 size selection protocol에 PCR 단계가 없다. HG03125 CLR은 historical freeze에만 남고 production에서는 1KGP Illumina로 대체됐다. |
| HGSVC3 | 60 | `pcr_free_protocol` | SMRTbell v3/Express 3 및 이전 HGSVC protocol에 PCR 단계가 없다. |
| GEUVADIS RNA | 449 pair | `pcr_amplified_confirmed` | TruSeq RNA Sample Prep v2에서 adapter-ligated cDNA를 10 PCR cycles로 enrichment한다. |

WGS 상태 합계는 `confirmed=2,918`, `protocol=62`, `inferred=232`다. RNA는 `pcr_amplified_confirmed=449`, WGS-only의 RNA 상태는 `not_applicable=2,763`이다.

## 실행 manifest 계약

`manifests/preprocess_manifest.tsv`와 각 `shard_manifest.tsv`에는 다음 열이 포함된다.

- `wgs_library_amplification`, `wgs_library_protocol`, `wgs_library_evidence`
- `rna_library_amplification`, `rna_library_protocol`, `rna_library_pcr_cycles`, `rna_library_evidence`

현재 production project는 `src/integrate_wgs_rna/library_provenance.py`의 canonical 값과 정확히 일치해야 한다. 누락 또는 불일치하면 입력 계약 검증에서 해당 sample만 실패하고 다음 sample로 진행한다. 사용자 추가 project는 모든 provenance 필드를 생략해 `unspecified`로 둘 수 있지만, 일부 필드만 채우는 것은 허용하지 않는다.

Stage 01 계약 JSON은 `mfm_wgs_input_contract_v2`, alignment QC는 `mfm_stage04_alignment_qc_v4`이며 성공·재사용·실패 모두 `library_provenance`를 기록한다. `stage04` schema 이름은 구형 내부 번호 호환이다.

## Duplicate 정책

- 1KGP Illumina aligned CRAM/BAM: 신뢰 가능한 upstream Picard MarkDuplicates provenance가 있으면 0x400 flag를 보존한다.
- raw 또는 upstream provenance가 없는 Illumina: Picard MarkDuplicates를 실행하되 `REMOVE_DUPLICATES=false`로 표시만 한다.
- HPRC/HGSVC PacBio: blanket Picard MarkDuplicates를 실행하지 않는다. source duplicate flag와 비율은 QC하며, 증폭/low-input 예외가 확인된 sample만 별도 opt-in 후보로 검토한다.
- GEUVADIS RNA: STAR 뒤 Picard MarkDuplicates를 실행해 표시만 하고, coverage/ASE read set에서 duplicate 제외 여부를 명시적으로 적용한다.

이 정책은 PCR 여부만으로 read를 물리 삭제하지 않는다. 정렬 결과와 read 수를 바꾸려면 별도 ablation과 승인이 필요하다.

## 생성·동기화·검증

```bash
conda run -n jj_mfm python scripts/build_portable_manifest.py
conda run -n jj_mfm python scripts/sync_shard_provenance.py
conda run -n jj_mfm python scripts/sync_shard_provenance.py --check-only
```

첫 명령은 main manifest와 summary를 만들고, 두 번째는 기존 shard 배정을 바꾸지 않은 채 provenance 열만 동기화한다. 세 번째는 37개 shard의 3,212명이 main manifest와 일치하는지 읽기 전용으로 확인한다.

## 근거

- 1KGP 30× WGS: https://pmc.ncbi.nlm.nih.gov/articles/PMC9439720/
- HPRC R1: https://www.nature.com/articles/s41586-023-05896-x
- HPRC R2 index: https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables/sequencing_data
- HGSVC2: https://pmc.ncbi.nlm.nih.gov/articles/PMC8026704/
- HGSVC3: https://www.nature.com/articles/s41586-025-09140-6
- GEUVADIS protocol: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1010958
- TOPMed RNA pipeline: https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md
