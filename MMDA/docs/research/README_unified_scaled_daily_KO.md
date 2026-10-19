# 2026-10-18 통합 논문 그림

원본 결과는 보존했다. 새 그림의 표시명은 아래첨자로 조판한
$\mathrm{DA}_{\mathrm{single}}$, $\mathrm{DA}_{\mathrm{multi}}$이다.
파일명과 CSV만 DA_single/DA_multi machine ID를 유지한다.
단일모델은 EnKF-NoPert, 다중모델은 기존 pooled MMDA-NoPert이다.

핵심 변경:
- Scaled SMAP daily로 filtered와 모든 SMAP skill/시계열/색상을 통일했다.
- 전체 성능 ACC는 모든 가용 연도의 circular 31-day climatology를 사용한다.
- Boxplot, triangle, 지도와 quadrant가 한 station 지표표를 공유한다.
- Boxplot의 평균 Delta ACC와 triangle의 평균 ACC 차이가 일치하는지 검증한다.
- All 356개, filtered 312개이다. 과거 raw-SMAP filtered는 317개였다.

## 그림 안내

| 파일군 | 내용 |
|---|---|
| fig00_landcover | Noah-MP 토지피복 |
| fig01_coverage_* | All/filtered network·site 분포, 일반 배경 및 토지피복 배경 |
| fig02_skill_* | ACC/raw ubRMSE/anomaly cRMSE boxplot |
| fig02_classification_* | 토지피복·평균 토양수분 분류별 dACC/dubRMSE |
| fig03_map_dACC*, fig03_map_dubRMSE* | OL/DA_single 대비 성능 변화 지도 |
| fig03_map_dCE*, dVE*, dTE* | CE/VE/TE 변화 지도, 원 단위 및 정규화 |
| fig04_quadrant_skill* | dACC-dubRMSE 사분면, scaled SMAP daily r 색상 |
| fig04_quadrant_CE_VE* | 오차 분해 사분면, DA_multi 자체 TE/nTE 색상 |
| fig_VE_CE_* | Triangle 분해 및 CE 아래 괄호 ACC |
| fig05~09 | 모델별 P, intra/inter, K/R, update, normalized innovation |
| fig10_* | 기존 3개 개선·3개 악화 사례의 daily/기후 anomaly/rolling ACC |
| fig11_monthly_* | 월별 전체 날짜 및 공통 scaled 관측일의 성능·P/K |

All과 filtered가 적용되는 모든 집단 그림은 두 버전으로 저장한다.
Grid P/K/innovation은 CONUS common-event 자료이며 station screening별로
전지구 grid를 다르게 자르지 않는다.

## 재현과 수치

- figures.ipynb: 모든 생성 코드.
- figures.executed.ipynb: 실제 실행 결과.
- METHODS_CAPTIONS_KO.md: 한글 계산법·수식·이전 ACC 불일치 원인.
- CAPTIONS_EN.md: 그림별 영어 caption.
- figure_manifest.csv: 전체 그림과 지표·표본·caption.
- delta_TE_color/: 기존 absolute TE 색상을 보존한 Delta TE 색상 추가판 8개.
  색상 부호는 사용자 요청대로 DA_multi minus baseline이므로 음수(파랑)가 개선이다.
- map_text_large.ipynb: 지도 내부 N/Mean/지점명 크기를 보완하는 지도 전용 코드.
- tables/station_metrics.csv: 유일한 전체 기간 성능표.
- tables/ACC_definition_audit.csv: 구 boxplot/triangle와 새 ACC 비교.
- tables/skill_summary.csv: All/filtered 평균.
- tables/paired_statistics.csv: boxplot 평균 개선·paired t-test.
- tables/station_year_month_metrics.csv: 월별 원 집계표.
- audit.json, verification.json, image_verification.csv: 실행·수치·시각 검증 기록.

처음부터 재계산하려면 새 kernel에서 figures.ipynb를 실행한다.
표기만 다시 그릴 때는 이 task에서 이미 계산한 월별 표를 다음처럼 재사용할 수 있다.
공유 cache나 다른 작업 폴더의 월별 성능표를 자동 재사용하는 옵션이 아니다.

~~~bash
MMDA_REUSE_LOCAL_MONTHLY=1 \
/home/gychoi/anaconda3/envs/kite/bin/python -m nbconvert \
  --to notebook --execute figures.ipynb --output figures.executed.ipynb \
  --ExecutePreprocessor.timeout=2400 --ExecutePreprocessor.allow_errors=False
~~~
