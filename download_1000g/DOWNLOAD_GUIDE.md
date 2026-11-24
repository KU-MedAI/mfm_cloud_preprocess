# Download guide

> **현재 적용 범위:** 이 가이드는 historical 3,202명 downloader에만 적용한다. 2026-09-16 PacBio 우선 3,212명 cohort는 `../share_package/DOWNLOAD_GUIDE.md`를 사용하며, PacBio run manifest가 동결되기 전에는 혼합 다운로드를 실행하지 않는다.

## 개요

이 문서는 `samples_3202.csv`와 `manifests/`만 사용해 1000 Genomes 3,202 sample 후보 데이터를 다운로드하는 방법을 정리한 가이드임.
모든 파일 경로는 이 package 폴더 기준 상대경로임.

## 1. 데이터 형식 요약

| 구분 | source/project | 파일 형식 | sample당 파일 수 | 바로 전처리 가능 여부 | 비고 |
|---|---|---|---:|---|---|
| Illumina short-read WGS | 1000 Genomes high-coverage 30x / ENA PRJEB31736 | `*.final.cram`, `*.final.cram.crai` | 보통 CRAM 1개 + CRAI 1개 | 가능 | GRCh38/hg38 aligned CRAM임. |
| ONT WGS aligned | 1KG_ONT_VIENNA | `*.hg38.cram`, `*.hg38.cram.crai` | 보통 CRAM 1개 + CRAI 1개 | 가능 | 현재 ONT 전처리에 가장 바로 쓰기 좋은 source임. |
| ONT WGS raw | HPRC ONT | BAM/uBAM | 여러 run 가능 | alignment 확인 필요 | aligned CRAM이 아니면 hg38/GRCh38 alignment 후 사용해야 함. |
| PacBio HiFi/CCS WGS | HPRC HiFi, HGSVC3 HiFi | BAM/uBAM | 여러 run 가능 | alignment 확인 필요 | sample별 모든 run을 받아야 coverage가 맞음. |
| RNA BAM | GEUVADIS TOPMed RNA-seq | STAR GRCh38 BAM | 보통 BAM 1개 + index | post-training pair 가능 | strict RNA-pair 기준 source임. |
| phased VCF | 1000 Genomes 2022 high-coverage phased panel | `*.vcf.gz`, `*.vcf.gz.tbi` | chr1-22 각각 | 공용 입력 | H1/H2 FASTA 생성과 haplotagging marker로 사용함. |

## 2. 기본 원칙

- `samples_3202.csv`는 sample-level 추천표임.
- 실제 다운로드는 `wgs_source_manifest`, `rna_source_manifest`에 적힌 `manifests/...` 파일을 기준으로 수행함.
- `wgs_primary_url`은 대표 URL임. Illumina/ONT Vienna처럼 sample당 alignment 파일이 1개인 경우에는 이 URL만으로 충분한 경우가 많음.
- HiFi/CCS, HPRC ONT는 multi-run일 수 있으므로 primary URL 하나만 받으면 안 됨. 반드시 manifest에서 같은 `sample`과 `source`의 모든 URL을 받아야 함.
- 대량 다운로드 전에는 `curl -I`, `wget --spider`, `aria2c --dry-run` 등으로 URL 상태를 확인하는 것이 안전함.

## 3. 준비

```bash
cd 3202_sample_share_package
mkdir -p downloads/wgs downloads/rna_bam downloads/vcf download_lists
```

대량 다운로드에는 `aria2c`를 권장함. 없는 경우 `wget -c -i`를 사용함.

```bash
aria2c -x 4 -s 4 -j 8 -c -i download_lists/urls.txt -d downloads/wgs
```

```bash
wget -c -i download_lists/urls.txt -P downloads/wgs
```

## 4. WGS 다운로드

### 4.1 Illumina 30x WGS CRAM

Illumina는 `manifests/illumina_30x_ena_final_cram_size_by_sample.tsv`를 사용함.
파일은 `*.final.cram`이며 index는 같은 URL에 `.crai`를 붙여 받음.

예시: 특정 sample 2개 다운로드 URL 생성.

```bash
python - <<'PY'
import csv
want = {"HG00096", "HG00097"}
with open("manifests/illumina_30x_ena_final_cram_size_by_sample.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["sample"] in want:
            url = row["path"]
            if not url.startswith(("http://", "https://", "ftp://")):
                url = "https://" + url
            print(url)
            print(url + ".crai")
PY
```

저장:

```bash
python - <<'PY' > download_lists/illumina_example.urls
import csv
want = {"HG00096", "HG00097"}
with open("manifests/illumina_30x_ena_final_cram_size_by_sample.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["sample"] in want:
            url = row["path"]
            if not url.startswith(("http://", "https://", "ftp://")):
                url = "https://" + url
            print(url)
            print(url + ".crai")
PY
aria2c -x 4 -s 4 -j 8 -c -i download_lists/illumina_example.urls -d downloads/wgs
```

### 4.2 ONT Vienna aligned CRAM

ONT Vienna는 `manifests/wgs_urls.tsv`에서 `source == ONT_Vienna`인 행을 사용함.
파일은 `*.hg38.cram`이며 index는 같은 URL에 `.crai`를 붙여 받음.

```bash
python - <<'PY' > download_lists/ont_vienna_example.urls
import csv
want = {"HG00096", "HG00100"}
with open("manifests/wgs_urls.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["sample"] in want and row["source"] == "ONT_Vienna":
            print(row["url"])
            print(row["url"] + ".crai")
PY
aria2c -x 4 -s 4 -j 8 -c -i download_lists/ont_vienna_example.urls -d downloads/wgs
```

### 4.3 HPRC ONT raw BAM/uBAM

HPRC ONT는 sample당 여러 run이 있을 수 있음.
`manifests/wgs_urls.tsv`에서 `source == HPRC_ONT`인 모든 행을 받아야 함.

```bash
python - <<'PY' > download_lists/hprc_ont_HG00099.urls
import csv
want = {"HG00099"}
with open("manifests/wgs_urls.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["sample"] in want and row["source"] == "HPRC_ONT":
            print(row["url"])
PY
aria2c -x 4 -s 4 -j 4 -c -i download_lists/hprc_ont_HG00099.urls -d downloads/wgs
```

주의:
- HPRC ONT는 `aligned_ref == -`인 raw/uBAM 성격의 파일이 많음.
- 현재 WGS 전처리에 바로 넣으려면 hg38/GRCh38 alignment 여부를 먼저 확인해야 함.

### 4.4 HPRC/HGSVC HiFi/CCS BAM/uBAM

HiFi/CCS도 sample당 multi-run인 경우가 많음.
sample/source별 모든 URL을 받아야 함.

```bash
python - <<'PY' > download_lists/hifi_HG00099.urls
import csv
want = {"HG00099"}
sources = {"HPRC_HiFi", "HGSVC_HiFi"}
with open("manifests/wgs_urls.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["sample"] in want and row["source"] in sources:
            print(row["url"])
PY
aria2c -x 4 -s 4 -j 4 -c -i download_lists/hifi_HG00099.urls -d downloads/wgs
```

주의:
- HPRC HiFi는 같은 sample에 여러 sequencing run이 존재할 수 있음.
- HGSVC HiFi도 여러 CCS/HiFi BAM으로 나뉠 수 있음.
- alignment가 되어 있지 않으면 minimap2/pbmm2 등으로 hg38/GRCh38 정렬 후 전처리해야 함.

## 5. RNA 다운로드

### 5.1 GEUVADIS processed RNA BAM

strict RNA-pair는 GEUVADIS processed RNA BAM 또는 ONT direct RNA BAM이 확인된 sample만 사용함.
GEUVADIS RNA BAM은 `manifests/rna_urls.tsv`에서 `source == Geuvadis_TOPMed`인 행을 사용함.

```bash
python - <<'PY' > download_lists/geuvadis_rna_bam_example.urls
import csv
want = {"HG00096", "HG00100"}
with open("manifests/rna_urls.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["sample"] in want and row["source"] == "Geuvadis_TOPMed":
            print(row["url"])
            if row["url"].endswith(".bam"):
                print(row["url"] + ".bai")
PY
aria2c -x 4 -s 4 -j 8 -c -i download_lists/geuvadis_rna_bam_example.urls -d downloads/rna_bam
```

### 5.2 recount3 GEUVADIS BigWig

recount3 BigWig은 이번 다운로드 대상에서 제외함.
원본 GEUVADIS RNA BAM이 아니라 재처리된 signal track이며, 현재 운영 기준은 RNA BAM pair임.
과거 BigWig audit manifest는 archive package에만 보존하고, v2 운영 폴더에는 포함하지 않음.

### 5.3 MAGE RNA FASTQ 제외

MAGE FASTQ는 RNA BAM이 아니므로 이번 운영 대상에서 제외함.

### 5.4 HPRC Kinnex RNA 제외

HPRC Kinnex/FLNC RNA는 RNA BAM 기준 운영 대상이 아니므로 이번 다운로드/전처리에서 제외함.

## 6. 공용 phased VCF 다운로드

공용 VCF는 2022 1kGP high-coverage phased SNV/INDEL/SV panel을 사용함.
chr1-22 각각 `vcf.gz`와 `vcf.gz.tbi`를 받음.

```bash
python - <<'PY' > download_lists/1kgp2022_phased_vcf.urls
base = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
for chrom in range(1, 23):
    url = f"{base}/1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    print(url)
    print(url + ".tbi")
PY
aria2c -x 4 -s 4 -j 8 -c -i download_lists/1kgp2022_phased_vcf.urls -d downloads/vcf
```

## 7. 그룹별 다운로드 URL 만들기

### 7.1 strict RNA-pair sample의 WGS만 받기

```bash
python - <<'PY' > download_lists/post_training_wgs.urls
import csv
with open("samples_3202.csv") as f:
    for row in csv.DictReader(f):
        if row["recommended_training_use"] == "post_training_wgs_rna_pair" and row["wgs_primary_url"]:
            print(row["wgs_primary_url"])
            idx = row["wgs_index_url_if_predictable"]
            if idx:
                print(idx)
PY
aria2c -x 4 -s 4 -j 8 -c -i download_lists/post_training_wgs.urls -d downloads/wgs
```

주의: 이 방식은 primary URL 중심임. HiFi/CCS 또는 HPRC ONT sample은 multi-run 누락이 생길 수 있으므로, 해당 sample은 `manifests/wgs_urls.tsv`에서 source별 전체 URL을 다시 뽑는 것이 안전함.

### 7.2 strict RNA-pair sample의 RNA BAM만 받기

```bash
python - <<'PY' > download_lists/post_training_rna_bam.urls
import csv
with open("samples_3202.csv") as f:
    for row in csv.DictReader(f):
        if row["recommended_training_use"] == "post_training_wgs_rna_pair":
            src = row["rna_download_source_available"]
            url = row["rna_primary_url"]
            if "Geuvadis_TOPMed" in src and url:
                print(url)
                if url.endswith(".bam"):
                    print(url + ".bai")
PY
aria2c -x 4 -s 4 -j 8 -c -i download_lists/post_training_rna_bam.urls -d downloads/rna_bam
```

### 7.3 WGS-only pre-training Illumina만 받기

```bash
python - <<'PY' > download_lists/pretraining_illumina_wgs.urls
import csv
with open("samples_3202.csv") as f:
    for row in csv.DictReader(f):
        if row["recommended_training_use"].startswith("pre_training") and row["assigned_wgs_platform"] == "Illumina_30x":
            if row["wgs_primary_url"]:
                print(row["wgs_primary_url"])
                if row["wgs_index_url_if_predictable"]:
                    print(row["wgs_index_url_if_predictable"])
PY
aria2c -x 4 -s 4 -j 8 -c -i download_lists/pretraining_illumina_wgs.urls -d downloads/wgs
```

### 7.4 WGS-only pre-training ONT Vienna만 받기

```bash
python - <<'PY' > download_lists/pretraining_ont_vienna_wgs.urls
import csv
with open("samples_3202.csv") as f:
    for row in csv.DictReader(f):
        if row["recommended_training_use"].startswith("pre_training") and row["wgs_download_source"] == "ONT_Vienna":
            if row["wgs_primary_url"]:
                print(row["wgs_primary_url"])
                if row["wgs_index_url_if_predictable"]:
                    print(row["wgs_index_url_if_predictable"])
PY
aria2c -x 4 -s 4 -j 8 -c -i download_lists/pretraining_ont_vienna_wgs.urls -d downloads/wgs
```

## 8. 다운로드 후 확인

### 8.1 파일 개수 확인

```bash
find downloads -type f | wc -l
```

### 8.2 CRAM/BAM index 확인

```bash
find downloads/wgs -name "*.cram" | while read -r f; do
  test -f "${f}.crai" || echo "[missing index] $f"
done
```

```bash
find downloads/rna_bam -name "*.bam" | while read -r f; do
  test -f "${f}.bai" || echo "[missing index] $f"
done
```

### 8.3 alignment 파일 quickcheck

```bash
find downloads/wgs downloads/rna_bam -name "*.cram" -o -name "*.bam" > download_lists/alignment_files.txt
samtools quickcheck -v $(cat download_lists/alignment_files.txt)
```

주의:
- `samtools quickcheck`는 기본적인 EOF/index 수준의 빠른 확인임.
- CRAM 내부 decode 오류나 reference mismatch는 별도 `samtools view` 또는 pipeline 단계에서 다시 확인해야 함.

## 9. 어떤 그룹을 받으면 되는가

| 목적 | 추천 필터 |
|---|---|
| pre-training WGS-only 전체 후보 | `recommended_training_use`가 `pre_training`으로 시작하는 sample |
| post-training WGS-RNA pair | `recommended_training_use == post_training_wgs_rna_pair` |
| 바로 사용 가능한 ONT WGS | `wgs_download_source == ONT_Vienna` |
| 바로 사용 가능한 Illumina WGS | `assigned_wgs_platform == Illumina_30x` |
| HiFi/CCS 후보 | `assigned_wgs_platform == HiFi_CCS` |
| RNA BAM pair | `has_rna_bam_strict == True` |

## 10. 주의사항

- RNA BAM이 실제 있는 strict pair는 449개임. RNA BAM이 없으면 모두 WGS-only로 통합함.
- recount3 BigWig은 이번 다운로드 대상에서 제외함. RNA-pair 판단은 RNA BAM 기준으로만 수행함.
- HiFi/CCS와 HPRC ONT는 multi-run이므로 source manifest에서 전체 URL을 받아야 함.
- HPRC raw/uBAM 계열은 전처리 전 alignment 상태를 반드시 확인해야 함.
- 모든 대량 다운로드는 sample subset으로 먼저 테스트한 뒤 확장하는 것이 안전함.


## 11. 운영 스크립트 빠른 실행 예시

### 11.1 기본 5 phase 생성

```bash
bash scripts/00_make_phase_assignments.sh
```

phase 수를 바꾸려면 다음처럼 실행함.

```bash
PHASES=10 bash scripts/00_make_phase_assignments.sh
```

### 11.2 phase 내부 shard 크기 조절

```bash
SHARD_GIB=750 bash scripts/01_make_phase_url_lists.sh 1
```

```bash
SHARD_GIB=2000 bash scripts/01_make_phase_url_lists.sh 1
```

### 11.3 source 선택 모드

aligned CRAM만 우선 다운로드함.

```bash
WGS_SOURCE_MODE=aligned_only RNA_SOURCE_MODE=strict_bam bash scripts/01_make_phase_url_lists.sh 1
```

sample별 대표 URL만 사용함.

```bash
WGS_SOURCE_MODE=primary RNA_SOURCE_MODE=primary bash scripts/01_make_phase_url_lists.sh 1
```

manifest의 모든 WGS run을 사용함.

```bash
WGS_SOURCE_MODE=all_manifest_runs bash scripts/01_make_phase_url_lists.sh 1
```

### 11.4 다운로드 dry-run

```bash
DRY_RUN=1 bash scripts/02_download_phase.sh 1 shard_001
```

### 11.5 실제 다운로드

```bash
bash scripts/02_download_phase.sh 1 shard_001
```

### 11.6 다운로드 없이 바로 전처리

이미 raw data가 준비되어 있으면 다운로드 없이 sample list를 만들고 전처리 command만 실행함.

```bash
MODE=pre PLATFORM=ONT bash scripts/05_preprocess_phase.sh 1 shard_001
```

```bash
RUN_PREPROCESS=1 \
PROCESS_CODE_DIR=../process_code \
PREPROCESS_CMD='SAMPLES="$SAMPLES_SPACE" ./run_cohort.sh igsr_ont_wgs_only' \
bash scripts/05_preprocess_phase.sh 1 shard_001
```
