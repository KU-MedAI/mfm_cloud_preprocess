# Readset → Pileup36 변경 코드

- 기존 readset과 새 pileup tensor의 변경 부분을 검토하기 위한 코드 원문·발췌본임.
- `code/01_pileup_원문`: 계산·36채널·VCF/좌표·저장·source 검증 원문임.
- `code/02_readset_발췌`: 비교 대상인 기존 readset 저장 부분임.
- `code/03_연결부_발췌`: 모드 선택·호출·reader·완료 검사 변경부임.
- `code/04_검사_원문`: 새 feature 수동 정답·저장·검증 테스트임.
- 출처·줄번호·복사/발췌 구분·SHA256은 `파일목록.csv`에 있음.

발췌본 단독 실행은 지원하지 않으며, 실제 실행은 전체 전달본의 진입점을 사용함.

## Feature 설명

- 신규 출력 순서: forward=3, reverse=4임. 36채널 수와 계산식은 유지하며 기존 순서도 기록된 feature_set으로 읽기·검증 지원함.

**[36개 feature의 계산식·분모·예시·해석](docs/pileup_36채널.md)** 한 문서에 모두 정리함.

개요 그림 4장과 상세 그림 18장으로 36채널 각각의 관측·분모·수치 예시를 설명함. 같은 지표의 공통/HP1/HP2는 나란히 비교함. 모든 숫자는 실제 계산 코드와 HTSlib 검증 경로로 확인한 합성 예시임.

![한 위치의 read를 세기](docs/figures/pileup36_easy/01_read_관측.png)

이전 흐름도·SVG 원본은 `figures/`에 보존함.

## 고정 CPU 활용 개선

- `code/05_CPU배분`: 고정 16CPU에서 HAP 작업을 겹쳐 실행하는 자원 계획과 검사 원문임.
- [적용 조건·스레드·RAM·검증 범위](docs/고정CPU_재정렬.md)에 정리함. 전체 실행 코드는 Cloud 저장소를 사용함.
