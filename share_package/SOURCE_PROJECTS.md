# 원천 프로젝트와 역할

| 프로젝트 | 이번 package의 역할 | 선택 인원 | 상태 |
|---|---|---:|---|
| HPRC R2 | PacBio HiFi read와 longcallD/DeepVariant 자료원 | 232 | cohort 선택 확정, file freeze 필요 |
| HGSVC3 | PacBio HiFi와 PAV SNV/indel 자료원 | 60 | cohort 선택 확정, file freeze 필요 |
| HGSVC2 | 선택 HiFi donor와 PAV Freeze4 자료원 | 1 | NA12878 선택; HG03125 CLR은 historical freeze만 보존 |
| HPRC R1 | HG03492 HiFi와 graph-derived VCF 자료원 | 1 | cohort 선택 확정, file freeze 필요 |
| 1KGP 30x | Illumina CRAM/CRAI와 phased cohort VCF | 2,918 | URL 확인, checksum/reference preflight 필요 |
| GEUVADIS TOPMed | processed RNA BAM/BAI | 449 | URL 확인, checksum preflight 필요 |
| HPRC R3 | membership provenance | 0 | production 제외 |

공식 root URL과 callset 설명은 `samples_3212.tsv`, `project_data_formats.tsv`, `manifests/public_vcf_sources.tsv`에 기록한다. 외부 전달 시에는 URL뿐 아니라 release/accession/checksum을 같이 고정해야 한다.
