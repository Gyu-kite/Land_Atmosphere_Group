# 2단계 실제 작업: 월별 자료동화 성능 비교 정합성 점검

작성일: 2026-10-19. 대상은 현재 연구에서 실제 사용하는 notebook이다.
이는 Skill의 새 대화 시험이 아니라 Skill을 추출하기 위한 실제 작업이다.

## 입력과 의도

실제 사용자 요청: “그럼 그다음엔 2단계를 해보면 되는건가 ? 2단계를 한번 해봐 그러면”.
후보 1을 실제로 수행하기 위해 에이전트가 정한 구체적인 작업은 다음과 같다.
아래 문장은 사용자 발언의 인용이 아니라 이번 작업의 명시적 범위이다.

> 연구 소개와 평가 기준을 읽고, monthly SMAP 비교 notebook에서 OL,
> DA_single, DA_multi, scaled SMAP의 기간·공통 날짜·표본·지표 정의·집계를
> 확인한다. 원본을 변경하지 말고, daily 시계열에서 가능한 수치를 재계산하여
> 일치·불일치·확인 불가를 근거와 함께 보고한다.

저장소 루트 README.md를 읽었다. README는 multi-source DA/LNETF 중심이므로
현재 MMDA의 지표·실험 조건은 대상 notebook과 그 METHODS 문서를 함께 읽었다.
2주차 제출 README의 정확한 사본인지는 확인되지 않았다.

대상 경로의 공통 prefix:

`notebooks/MMDA_paper_rebuild_20261004/main/new_20261018/unified_scaled_daily/`

- 대상 notebook: `monthly_sm_climatology/smap_skill/figures.ipynb`
- 대응 Python: `monthly_sm_climatology/smap_skill/figures.py`
- 실행 기록: `monthly_sm_climatology/smap_skill/figures.executed.ipynb`
- 그림의 값: `monthly_sm_climatology/smap_skill/monthly_same_date_skill.csv`
- 상위 계산: `figures.py`의 `clim31`, `components`, monthly 계산 부분
- 상위 수치표: `tables/station_metrics.csv`, `tables/station_year_month_metrics.csv`
- 분류: `monthly_sm_climatology/station_membership.csv`

## 실제 수행 절차

1. notebook에서 읽는 CSV만으로 결론을 내리지 않고 상위 계산 코드와 cache
   생성 코드까지 추적했다.
2. scaled SMAP의 출처가 `mmda_longterm_qcv4_paper.py`의 EnKF DAOBS
   두 번째 record와 UTC daily resampling이라는 것을 코드에서 확인했다.
3. daily OL/DA_single/DA_multi/ISMN cache와 scaled SMAP cache를 읽었다.
4. 같은 물리값 날짜를 먼저 선택하고, 별도 구현으로 모든 가용 연도의
   circular 31-day climatology를 계산했다. 단순 time-axis rolling이 아니다.
5. 관측일 공통 anomaly에서 월당 최소 8일 조건을 적용하고, ACC와 absolute
   daily error의 ubRMSE를 직접 재계산했다.
6. filtering을 absolute daily Pearson r 차이 > -0.2와 최소 60일 조건으로
   재확인했다. 이는 월별 anomaly ACC와 다른 정의이다.
7. 지점/달마다 연도 평균 후 지점 평균하는 집계를 다시 계산했다.
8. 입력 파일 SHA-256을 실행 전후 비교해 원본 비변경 여부를 검사했다.

점검 스크립트는 대상 plotting helper 함수를 그대로 호출하지 않는다.
다만 같은 daily cache를 사용하므로 원자료까지 독립된 검증은 아니다.

## 점검표

| 항목 | 판정 및 확인 범위 | 근거 |
|---|---|---|
| SMAP 종류 | 코드상 scaled field 사용 확인 | 상위 figures.py의 `['scaled']`; cache builder의 DAOBS record 처리 |
| 비교 날짜 | 공통 관측일에서 계산하는 구조 확인; 수치 재검산 결과는 audit.json | 상위 monthly의 `dropna()`와 재계산 date_support.csv |
| ACC 정의 | all-year circular 31-day anomaly correlation | 상위 clim31 및 별도 anomaly 구현 |
| ubRMSE 정의 | 동일 날짜 absolute daily error의 population std, ddof=0 | 상위 monthly 계산 및 재계산 |
| Filtering | absolute daily Pearson 차이 > -0.2, 최소 60일 | station_metrics와 daily screening 재계산 |
| Monthly 집계 | 연도 → 지점/달 → 지점 동일 가중 평균 | 대상 figures.py:38–43, group_output 재계산 |
| 차이 부호 | ACC는 j−OL, ubRMSE는 OL−j로 양수가 개선 | 대상 figures.py:42–43 |
| 지점 분류 | 고정 OL mean SM tertile, All/filtered 각각 경계 | membership와 station_metrics |
| SNOTEL | 최종 station table에 없음 | network 검사 |
| 실제 L2/DAOBS와 cache 일치 | 확인 불가: 이번 실행은 원자료 전체를 다시 읽지 않음 | cache builder 코드만 추적 |
| MMDA 고유 scaled 관측 | 이 그림에서는 사용하지 않음 | 공통 SMAP은 DA_single에서 유래 |
| 동화 지원 metadata 자체의 정확성 | 이번 범위에서 독립 재검증하지 않음 | 상위 eligibility metadata 재사용 |

수치 일치 여부의 최종 결과와 최대 차이는 `audit.json`에 기록한다.
실행이 실패하면 그 기록을 완료 결과로 대체해서는 안 된다.

### 이번 실행의 최종 수치

점검 스크립트는 정상 종료(exit 0)했다. 356개 지점의 관측일 월별 조합
20,401개를 재계산하고, 최종 그림 집계 96행과 비교했다.

| 검사 | 최대 절대 차이 | 결과 |
|---|---:|---|
| 월별 ACC | 8.07 × 10^-14 | 일치 |
| 월별 ubRMSE | 9.97 × 10^-17 | 일치 |
| TE − (CE + VE) | 1.39 × 10^-17 | 일치 |
| Screening correlation 차이 | 9.71 × 10^-17 | 일치 |
| 최종 group 평균 지표 | 3.33 × 10^-16 | 일치 |

지점·연도·월별 N, 월별 집계 N, 고정 SM class도 검사에 통과했다.
검사한 입력의 실행 전후 SHA-256이 모두 같아 원본 비변경을 확인했다.
따라서 이번 확인 범위에서는 수정할 수치 불일치가 발견되지 않았다.

## 해석 시 유의할 사항

### 월별 표본은 고정되어 있지 않다

기본 등록 표본은 All 356, filtered 312개이다. 하지만 월별 공통 관측일 기준을
만족하는 실제 지점 수는 All 173–349, filtered 153–306개이다.
이 수의 변화 자체는 구현 오류가 아니지만, 계절별 차이를 동일 지점의
순수한 계절 효과로만 해석해서는 안 된다. 고정 지점 집단의 계절 효과가
목적이라면 별도 common-station sensitivity 분석이 필요하다.

### SMAP skill의 의미를 구분한다

SMAP은 model-space scaled daily 관측이며 satellite native retrieval skill이 아니다.
모든 선의 날짜는 맞추지만 SMAP 자체를 각 실험의 고유 scaled 자료로
각각 평가한 그림은 아니다. 현재 그림은 공통 DA_single-scaled SMAP과
실험 결과를 비교하는 그림이다.

### 검증 범위를 확대해서 말하지 않는다

daily cache에서 월별 수치가 재현되어도 모델 run, 관측 reader, 원자료 QC,
개별 동화 시각의 배치가 모두 정상이라는 뜻은 아니다. climatology에는
평가 연도가 포함되므로 out-of-sample 예측 검증으로 표현하지 않는다.

## 수정 요청과 초안 추출 기록

이번 단계의 최초 결과 이후 사용자의 새 수정 요청은 아직 없다.
과거의 다음 요청을 이미 알려진 요구사항으로 적용했으며, 이번 단계에서
새로 받은 피드백인 것처럼 기록하지 않는다.

- “SMAP을 scaled로 바꿔야 한다”: raw/scaled 출처를 코드에서 확인.
- “climatological 31-day mean이다”: time-axis rolling과 구분.
- “CSV에 의존하지 말고 독립적으로 계산”: daily cache에서 수치 재계산.
- “SMAP과 같은 날짜로 비교”: 다섯 series 공통 유효 날짜를 사용.

이번 실제 수행에서 추출한 Skill 초안은 `SKILL_DRAFT.md`이다.
정식 `.agents/skills/` 설치본이 아니며, 사용자 직접 검토를 아직 받지 않았다.
특히 추출된 교훈은 “코드상 정상”과 “수치 재현”과 “원자료 검증”을 구분하고,
표본 변화를 반드시 보고하며, 모르는 조건을 추정하지 않는 것이다.

## 재현 명령 및 산출물

저장소 루트에서:

```bash
/home/gychoi/anaconda3/envs/kite/bin/python \
  docs/assignement/step2_audit_20261019/audit_comparison.py
```

- `audit.json`: 파일 경로·SHA-256, 재계산 오차, 검증 한계
- `recalculated_station_year_month.csv`: 독립 재계산 월별 지표
- `date_support.csv`: 지점별 유효 날짜 수·시작·종료
- `sample_counts.csv`: All/filtered 및 SM class별 실제 월별 표본
- `SKILL_DRAFT.md`: 이번 절차를 좁게 추출한 미검토 Skill 초안

기존 모델·관측·cache·notebook·그림은 수정하지 않았다.
