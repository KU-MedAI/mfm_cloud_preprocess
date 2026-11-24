# v2 전처리 파이프라인 교육자료

## 제공 파일

- `전처리_파이프라인_교육자료.pptx`: 16:9, 38장
- `전처리_파이프라인_교육자료.pdf`: LibreOffice 렌더링 고정본
- `교육자료_슬라이드원고.md`: 슬라이드별 설명과 근거
- `교육자료_근거맵.tsv`: 슬라이드–문서/코드 매핑
- `preview/contact_sheet.png`: 전체 슬라이드 미리보기
- `preview/slides/`: 38장 개별 미리보기
- `build_training_ppt.py`: PPTX 재생성 코드
- `requirements.txt`: 생성 의존성
- `assets/`: 2026-09-28 도식 원본 보존용. 현재 PPT는 오래된 번호가 다시 들어오지 않도록 flow를 코드에서 직접 그림

## 대상과 사용법

현재 v2를 처음 인수받는 분석가·모델 설계자·운영 담당자를 위한 약 50~60분 교육 자료다. 발표자는 먼저 PPT의 발표자 노트를 읽고, 세부 근거가 필요할 때 `교육자료_슬라이드원고.md`와 `교육자료_근거맵.tsv`를 사용한다.

자료는 2026-10-02의 v2 코드·3,212명 manifest·Markdown을 기준으로 한다. 실행 단계는 WGS 01–12, RNA 13–26, cohort catalog 27을 기준으로 하며, 과거 Stage 22/70 명칭은 호환성 설명으로만 병기한다.

## 재생성

`jj_mfm` 환경의 고정 의존성을 사용한다.

```bash
conda activate jj_mfm
python build_training_ppt.py

mkdir -p render
libreoffice --headless --convert-to pdf --outdir render \
  전처리_파이프라인_교육자료.pptx
```

## 검증된 상태

- PPTX ZIP 무결성: PASS
- PPTX/PDF 페이지: 38장
- 별도 슬라이드 원고와 근거맵: 생성 완료
- 로컬 문서·코드 근거 파일명: 누락 0건
- 슬라이드 경계를 벗어난 shape: 0건
- PDF 16:9 렌더링 및 전체 contact sheet 육안 검토: PASS

## 업데이트 원칙

단계 번호, 필터 preset, 표본 수, 병렬도 또는 manifest schema가 바뀌면 PPT만 직접 수정하지 말고 `build_training_ppt.py`와 근거맵을 함께 수정한 뒤 PPTX·PDF를 다시 생성한다.
