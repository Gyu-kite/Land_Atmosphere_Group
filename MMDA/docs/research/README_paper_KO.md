# MMDA 논문 통합 분석

## 현재 투고용 구성: 2026-10-05

**최종 그림 선택과 순서는 `submission_filtered_20261005/README_KO.md`를 따른다.**
Filtered를 main 결과로, All과 기존 Fig11·12를 supplementary로 옮겨 복사했다.
원래 main/appendix/timeseries는 재현 및 보존용이며, 아래 과거 Fig 번호는 원본 번호다.
새 구성은 10개 figure 번호(분리 panel 포함 12개 파일), 보충자료 32개 파일이다.
`submission_filtered_20261005/figure_manifest.csv`에 새 번호와 원본 파일의 대응이 있다.
`organize_submission.py`로 삭제 없이 다시 구성할 수 있다. Plot 코드는 기존 통합 notebook에 유지한다.

## 실행 및 편집

`MMDA_paper_all_figures.ipynb`에 성능 평가부터 intra/inter 메커니즘까지 모든 계산 및 그림 코드가 들어 있다. 별도 메커니즘 노트북 실행이나 `%run`은 필요 없다. 공통 daily/event 입력 cache는 기존 장기 분석 cache를 사용한다. 모델 output을 변경하거나 LIS를 재구동하지 않는다.

- Notebook: `MMDA_paper_all_figures.ipynb`
- 대응 source: `build_figures.py`
- 전체 수치·그림 점검: `verify_outputs.py`, `verify_mechanism.py`
- 편집 전 원본: `archive/pre_integration_20261004/`

노트북을 직접 수정한 후에는 `make_notebook.py`로 덮어쓰지 않는다. 이 도구는 source에서 노트북을 다시 만드는 용도다.

## 논문 흐름

| 질문 | 그림 |
|---|---|
| 어디에서 어떤 표본을 비교했는가? | Fig01 coverage |
| EnKF-NoPert 및 OL보다 개선됐는가? | Fig02-04 skill, 지도, 사분면 |
| 어떤 모델의 spread와 불일치가 P를 구성했는가? | Fig05-06 모델별 P, intra/inter |
| 실제 관측 반응은 어떻게 달라졌는가? | Fig07-09 gain/R, update, innovation |
| 평균 gain만 같게 해도 차이가 남는가? | Fig11-12 one-step 대조, 독립 cycling 실험 아님 |
| 현재 update 전/후 오차 차이는 어디서 생겼는가? | Fig13 exact centered-error accounting |
| Inter 추가 update는 항상 도움이 되는가? | Fig14 이득, 잘못된 방향, overshoot |
| 계절별 크기와 이득은 일치하는가? | Fig15 동일 station 집합의 월별 비교 |
| 개선/악화 사례에서 실제 변동은 어떤가? | Fig10 기본 시계열 + Fig16 같은 사례의 update 상세판 |
| 대조군 및 진단 한계는? | FigA01 screening, A02 EnKF-Pert, A03 event 구성, A04 후속 cycling |

Fig10과 Fig16은 서로 다른 station을 추가한 것이 아니라 동일 사례를 skill 관점과 update 관점에서 연결한 것이다. 최종 투고 시 지면에 따라 Fig16 일부를 보충자료로 배치할 수 있다.

## 결과 읽기

`RESULTS_KO.md`는 전체 논문 해석, `MECHANISM_RESULTS_KO.md`는 세부 수치와 실패 사례다. `CAPTIONS_METHODS_KO.md` 및 `MECHANISM_METHODS_KO.md`에 지표·수식을 설명한다. 전체 그림별 한글 캡션은 `FIGURE_CAPTIONS_KO.md` 하나로 관리한다.

추가 표와 재현용 event 자료는 `mechanism/tables/`, `mechanism/cache/`에 있고 그림은 기존 `main/`, `appendix/`, `timeseries/`에 함께 저장한다. 별도의 `mechanism/figures`는 만들지 않는다.

## 해석의 경계

Inter가 더한 update의 평균적인 유익함과 inter가 단순 inflation보다 고유하게 우수하다는 주장은 다르다. 후자는 독립 cycling 대조 없이 확정하지 않는다. 순간 update에는 daily ISMN을 proxy로 사용한다. Event cMSE 분해를 daily ACC의 인과 기여율로 바꾸지 않는다. 이후 시계열에는 추가 동화가 포함될 수 있다.
