# MMDA: Multi-Model Soil Moisture Data Assimilation

LISF 기반 토양수분 Multi-Model Data Assimilation(MMDA)의 연구 안내와
반복 작업 Skill 과제를 정리한 문서 패키지이다. 전체 LISF 소스나 관측·모델
대용량 자료를 포함하는 배포판이 아니다.

## 연구 개요

Noah-MP 5.0, Noah 3.9, CLSM의 ensemble을 reference space에서 결합하고,
pooled background variance를 사용하는 자료동화를 단일모델 EnKF와 비교한다.
ISMN QC-v4로 daily 토양수분을 평가하며 SNOTEL은 제외한다.
주요 지표는 climatological 31-day anomaly correlation과 ubRMSE이다.
Pooled variance의 within-model 및 between-model 성분을 분리하되,
모델 간 불일치를 검증된 실제 오차분산이나 고유한 인과 기여로 단정하지 않는다.

## Skill 과제

- [정식 Skill](.agents/skills/audit-da-comparison/SKILL.md)
- [연구 배경과 평가 기준](.agents/skills/audit-da-comparison/references/research.md)
- [5주차 보고서](docs/assignement/05주차_과제_보고서.md)
- [제출 브랜치와 수정 이력](SUBMISSION_GIT.md)

Skill은 자료동화 평가 notebook의 날짜·표본·관측 scaling·지표 정의·집계를
점검한다. 모델 구동이나 알고리즘 변경은 수행하지 않는다.
실제 입력 notebook과 물리값 cache는 원래 연구 환경에 있으며 이 패키지에
포함하지 않는다. 제출한 감사 코드·수치표는 수행 근거이지 전체 재실행 환경은 아니다.
저장소를 내려받은 뒤 MMDA/를 작업 디렉토리로 열면 그 안의
.agents/skills/가 프로젝트 Skill 경로가 된다.

## 복사해 보존한 연구 안내

- [통합 논문 분석 안내](docs/research/README_paper_KO.md): 2026-10-04 문서 계열.
- [최신 통합 평가 안내](docs/research/README_unified_scaled_daily_KO.md): 2026-10-18 개정.
- [Patch testbed 안내](docs/research/README_patch_testbed.md).
- [Patch 실행 흐름](docs/research/README_testbed_workflow.md).

복사 문서 안의 상대 경로는 원래 연구 디렉토리 기준이다. 연결된 notebook,
cache, 모델 output 모두를 이 패키지로 옮긴 것은 아니다. 최신 평가 정의는
2026-10-18 개정 문서와 Skill reference를 먼저 확인한다.
기존 README와 원본 연구 파일은 삭제하거나 이동하지 않고 복사했다.
