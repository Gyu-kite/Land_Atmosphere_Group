---
name: audit-da-comparison
description: Audit soil-moisture data-assimilation evaluation notebooks for comparable dates, samples, observation scaling, metric definitions, and aggregation. Use when checking whether experiment skill comparisons are consistent; not for model runs, figure styling, or algorithm design.
---

# 자료동화 성능 비교 정합성 점검

## 입력과 범위

평가 notebook 또는 Python 코드, 실제 자료 경로, 의도한 비교 기준을 받는다.
먼저 [연구 배경과 기준](references/research.md)을 읽는다. 사용자가 별도 연구
소개를 제공하면 함께 읽고, 현재 요청과 과거 기본값을 구분한다.
필수 경로나 비교 기준이 없고 코드에서도 확정할 수 없다면 추정하지 말고 질문한다.
확인 가능한 부분만 먼저 점검하고 미확인 항목을 남긴다.

## 절차

1. notebook의 실제 입력과 상위 계산 코드를 추적한다. 표시 이름만으로
   실험, 관측 단위, raw/scaled 여부를 판단하지 않는다.
2. 기간, 공통 날짜, 지점 및 screening 조건을 확인한다. 각 지표의 실제
   유효 표본 수와 지점 구성을 확인하고, 표본이 다르면 계절 효과와 구분해 보고한다.
3. ACC, ubRMSE, climatology, 정규화와 차이 부호를 구현에서 확인한다.
   코드와 명시된 기준이 다르면 불일치로 보고하고 임의로 기준을 바꾸지 않는다.
4. 가능한 범위에서 원시 시계열 또는 가장 가까운 상위 자료로 수치를 다시
   계산한다. 동일 함수 재호출과 독립적인 수치 재계산을 구분한다.
5. 동일 표본의 평균 차이, 성분 합, 집계값을 검산한다. 재계산 범위와 최대 차이를
   기록하고, 수치 허용 오차는 단위·자료 정밀도에 맞게 정한다.
6. 원본을 보존하고 요청된 출력 또는 별도 작업 경로에 근거를 저장한다.
   수정 요청이 없으면 notebook·원자료·공유 cache를 변경하지 않는다.

## 출력 순서

1. 대상과 비교 기준: 실제 파일, 비교 실험, 지표 및 평가 기간.
2. `일치 / 불일치 / 확인 불가` 점검표: 실제 값, 코드·자료 위치, 확인 범위.
3. 발견 사항과 영향: 오류와 단순한 해석 한계를 구분.
4. 필요한 수정 또는 질문: 수치 수정이 불필요하면 그렇다고 명시.
5. 재현 정보: 실행 명령, 결과 파일, 재계산 범위와 남은 한계.

## 해석과 실행 경계

- cache 수준의 재현을 원자료 QC나 모델 구동 전체의 검증으로 확대하지 않는다.
- 읽어보지 않은 자료를 정상이라고 판정하거나 실행하지 않은 시험을 완료라고 쓰지 않는다.
- 관측 skill과 DA 개선의 연관성을 인과 효과로 쓰지 않는다.
- 모델 실행, 컴파일, 알고리즘 변경, 그림 스타일 수정은 이 Skill의 범위가 아니다.
