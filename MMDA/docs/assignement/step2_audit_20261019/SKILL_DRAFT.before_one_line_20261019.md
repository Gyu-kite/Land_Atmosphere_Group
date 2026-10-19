---
name: audit-da-comparison
description: Audit soil-moisture data-assimilation evaluation notebooks for comparable dates, samples, observation scaling, metric definitions, and aggregation. Use when checking whether reported experiment comparisons are consistent; not for model runs, figure styling, or algorithm design.
---

# 자료동화 비교 정합성 점검

평가 notebook, 실제 입력 경로, 의도한 비교 기준을 받아 정합성을 점검한다.
기준이 없거나 자료가 모호하면 추정하지 말고 질문한다. 사용자 요청 없이
실험 조건이나 원본 notebook을 수정하지 않는다.

1. notebook의 실제 입력과 상위 계산 코드를 추적한다. 표시 이름만으로
   실험, 관측 단위, raw/scaled 여부를 판단하지 않는다.
2. 기간, 공통 날짜, 지점 및 screening 조건을 확인한다. 각 지표의 실제
   유효 표본 수를 확인하고, 월별 표본 변화를 누락하지 않는다.
3. ACC, ubRMSE, climatology, 정규화와 차이 부호를 구현에서 확인한다.
   연구별 기본값은 references/research.md가 있을 때 참조하되 명시적 요청을 우선한다.
4. 가능한 범위에서 원시 시계열 또는 가장 가까운 상위 자료로 수치를 다시
   계산한다. 동일 함수 재호출과 독립적인 수치 재계산을 구분한다.
5. 동일 표본의 차이·평균·성분 합을 검사한다. 원자료까지 확인하지 못했다면
   그 범위를 밝히고 전체 파이프라인이 정상이라고 확대 해석하지 않는다.
6. 원본을 보존하고 별도 경로에 근거를 저장한다. 결과를 다음 순서로 보고한다.

## 출력 순서

- 대상과 비교 기준
- `일치 / 불일치 / 확인 불가` 점검표: 실제 값, 코드·자료 위치, 확인 범위
- 발견 사항과 결과에 미칠 영향
- 필요한 수정 또는 사용자에게 물어야 할 사항
- 실행 명령, 결과 파일, 재계산 범위와 남은 한계

오류가 없으면 없다고 쓰되, cache 수준 검증과 원자료 검증을 구분한다.
성능 개선을 만들기 위해 기준을 바꾸거나 관측 skill과 개선의 상관을 인과로 쓰지 않는다.
