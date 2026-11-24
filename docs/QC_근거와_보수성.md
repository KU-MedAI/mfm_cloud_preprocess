# QC 기준의 근거와 보수성

이 문서는 실행 기준의 출처를 명확히 하기 위한 전달본 요약이다. 전체 stage별 기준은 `단계별_입출력_상세.md`를 따른다. 여기서 **문헌의 원칙**과 **이 pipeline의 운영값**을 구분한다.

| 구분 | 대표 기준 | 처리 |
|---|---|---|
| 무결성 guardrail | 파일/index, `samtools quickcheck`, sample ID, FASTA/VCF/CRAM 22-autosome LN·M5·REF, chain, trace, phased tag GT, producer→consumer SHA-256, final package | 숫자를 완화하지 않는다. 불일치는 해당 sample만 중단하고 다음 sample으로 넘어간다. |
| 공식 workflow | paired FASTQ→STAR→Picard MarkDuplicates, `sjdbOverhang=read length-1`, unique STAR alignment의 MAPQ 255 | TOPMed/GTEx 구현을 따른다. 개인 HAP FASTA/lifted GTF만 본 pipeline의 목적에 맞게 다르다. |
| calibration 대상 | WGS HAP mapped fraction, RNA paired fraction, MAPQ 30, window N/blacklist, CNV rd 0.5–1.5 | 문헌의 방향은 있으나 보편적 동일 숫자는 아니다. 기본은 관측값과 기존 0.90 참고선 대비 상태만 QC에 기록하고 sample을 중단하지 않는다. cohort distribution·retained read/depth·allele balance로 재검토한다. |
| 모델 heuristic | ASE split-half, informative fraction, eligible gene 수, Pearson, overdispersion fit | 전처리 성공과 분리된 `sample_training_ready` label이다. 원본/최종 WGS·RNA package cleanup integrity와 혼동하지 않는다. |

## 비율 QC 기본 정책

`0.90`은 논문이나 TOPMed의 universal QC cutoff가 아니다. 따라서 기본 `QC_RATIO_POLICY=record_only`에서는 WGS HAP mapped fraction과 Illumina paired fraction, RNA paired flag fraction과 FASTQ pair 회수율을 수치로만 남긴다. 0.90 미만이어도 이 비율만으로 sample을 중단하거나 read/window를 제거하지 않는다.

과거 hard-fail 동작을 정확히 재현해야 할 때만 `QC_RATIO_POLICY=strict`를 명시한다. strict에서는 기존 `WGS_MIN_MAPPED_FRAC=0.90`, `RNA_MIN_PAIRED_FRACTION=0.90`을 중단 기준으로 쓴다.

- **Stage 02 input:** `paired primary records / all primary records`. paired SAM flag 비율이며 proper-pair completeness 자체는 아니다.
- **Stage 13 FASTQ:** `2×R1/R2 pairs / (2×pairs + singleton discard + unpaired discard)`. R1=R2 및 pair 수>0도 동시에 요구한다.

두 단계의 분모가 다르므로 두 결과와 `below_reference`, `qc_action`을 각각 QC JSON에 남긴다. 다만 R1/R2 개수 불일치, pair 0개, 서열·품질 손상 같은 구조적 문제는 비율 QC와 다르므로 계속 차단한다. single-end RNA는 현재 pair/fragment 기반 ASE 모델과 맞지 않으므로 별도 branch 설계 없이는 수용하지 않는다.

## 새 무결성 검사는 왜 hard-fail인가

- Illumina tag VCF의 모든 GT가 `0|1` 또는 `1|0`인지 확인하는 것은 allele 비율 cutoff가 아니다. unphased/partial GT는 HAP1/HAP2 방향 자체를 정의하지 못하므로 해당 sample을 중단한다. `fast_validate`도 이 계약은 생략하지 않는다.
- count·CNV·coverage SHA-256은 좋은 생물학적 신호의 비율을 판정하지 않는다. Stage 21/22가 만든 바로 그 내용을 Stage 23/25가 읽었는지, 검증 뒤 파일이 바뀌지 않았는지를 확인한다.
- execution provenance와 stage log 검사는 WGS 01–12 및 RNA 13–24가 성공 또는 명시적으로 skip되었음을 확인한다. 임의 threshold 미달을 새 실패 기준으로 추가하지 않는다.
- 따라서 이 항목들은 `record_only` 비율 QC와 달리 완화 대상이 아니다. 실패해도 batch 전체가 아니라 해당 sample만 차단한다.

## 출처

- TOPMed/GTEx RNA pipeline: paired 2×101 input, STAR, Picard MarkDuplicates, FASTA index 및 index parameter. <https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md>
- ENCODE bulk RNA standards: paired/single-end 구분과 paired library의 read-pair quantification. <https://www.encodeproject.org/data-standards/encode4-bulk-rna/>
- 1KGP 3,202 high-coverage WGS: BWA-MEM, Picard MarkDuplicates, BQSR functional-equivalence workflow. <https://pmc.ncbi.nlm.nih.gov/articles/PMC9439720/>
- ASEP: total count 20 및 minor allele count 5 filter, 5% FDR. <https://pmc.ncbi.nlm.nih.gov/articles/PMC7241832/>
- SAM specification: MAPQ의 Phred-scale semantics; aligner별 calibration이 필요하다. <https://samtools.github.io/hts-specs/SAMv1.pdf>
