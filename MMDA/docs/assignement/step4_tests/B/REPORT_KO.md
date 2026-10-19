# 월별 SMAP·자료동화 성능 비교 점검

2026-10-19에 README.md, 설치된 audit-da-comparison 스킬과 연구 기준을 읽고 수행했다. 루트 README는 LNETF/multi-source 개요이며, 이 notebook의 구체적 비교 기준은 상위 계산 코드와 METHODS_CAPTIONS_KO.md에서 확인했다.

## 대상과 기준

대상: notebooks/MMDA_paper_rebuild_20261004/main/new_20261018/unified_scaled_daily/monthly_sm_climatology/smap_skill/figures.ipynb.

입력 추적: 대상 notebook → tables/station_year_month_metrics.csv 및 monthly_sm_climatology/station_membership.csv → 상위 unified_scaled_daily/figures.py → notebooks/mmda_longterm_qcv4_paper.py의 daily cache 생성 코드. 실제 파일명과 SHA-256은 audit.json에 있다.

실험 경로는 /land4/LSM_OUTPUT/MMDA 아래 OL=MP50_OL_init20141231_through202512_direct18c, DA_single=MP50_EnKF30_nopert_2014_2025, DA_multi=MP50_39_CLSM_Super_Pooled_anomcdf_2018_2025이다. 디렉토리명의 연도로 평가 기간을 추정하지 않았다. 모델 변수는 SoilMoist_tavg이다. 이 감사는 실험 설정 파일과 ensemble 구성 자체를 재검증하지 않았다.

## 점검표

| 항목 | 판정 | 실제 값·확인 범위 및 근거 |
|---|---|---|
| 평가 날짜 | 일치 | 상위 figures.py:530–546에서 다섯 anomaly series의 dropna 후 월별 계산. cache builder의 평가 축은 2015-01-02–2025-12-31. 실제 날짜 범위는 supplement.json, 지점별 범위는 date_support.csv |
| 방법 간 유효 표본 | 일치 | 월당 최소 8 공통일, 비영 anomaly 분산. 모든 관측일 station-year-month의 N을 daily cache에서 재확인 |
| 월간 고정 표본 | 불일치 | 각 월은 같은 지점 집합이 아니다. 이는 구현 오류가 아니라 계절 해석의 한계. monthly_site_composition.csv에 실제 key 및 전월 대비 교체를 기록 |
| 관측 종류 | 일치 | 평가 기준은 qc_v4 ISMN, SNOTEL 제외. SMAP은 EnKF DAOBS 두 번째 record의 scaled daily mean. native retrieval이나 MMDA 고유 관측을 혼합하지 않음 |
| 각 run의 실제 동화 센서·QC·update 표본 | 확인 불가 | reader 및 run 설정, 원본 DAOBS 전량 재검증은 수행하지 않음. 위 관측 종류 판정은 평가 입력에 한정 |
| ACC | 일치 | 상위 figures.py:49–62, 64–87. 전 가용 연도 calendar-day 평균 후 순환 ±15일, 최소 15 유효 bin. Feb29 학습 제외, Feb28/Mar1 보간. anomaly Pearson 상관 |
| climatology의 학습 표본 | 방법 간 다름, 문서와 일치 | 모델/ISMN은 네 series 공통 일자, SMAP은 가용 관측일. 평가 날짜는 같지만 climatology 학습 날짜까지 같다는 뜻은 아님. METHODS_CAPTIONS_KO.md:110–111 |
| ubRMSE | 일치 | 상위 figures.py:546에서 공통 날짜의 absolute daily error 표준편차, ddof=0. 절댓값을 취한다는 뜻이 아니라 anomaly를 빼기 전 물리값의 오차임 |
| screening | 일치 | ≥60 paired days, 비영 분산, r(ISMN,scaled SMAP)−r(ISMN,OL)>−0.2. anomaly ACC가 아닌 일평균 Pearson r |
| 공간·지점 기준 | 일치 | assim_supported & in_conus, ≥365 common anomaly days. membership은 추가 lc_code>0; All은 모든 원시 ISMN 지점이라는 뜻이 아님 |
| 월별 집계 | 일치 | 대상 figures.py:38–43: 지점/달 내 연도 평균 후 지점 동일 가중. 전체 평균은 세 class 평균의 단순 평균이 아님 |
| SM class | 일치 | All/filtered 각각 OL_mean_SM tercile로 고정. 두 subset의 경계는 같다고 가정할 수 없음 |
| 개선 부호 | 일치 | ΔACC=j−OL, ΔubRMSE=OL−j. 양수 개선. 최종 CSV 차이 열도 재검산 |
| 정규화 | 해당 없음 | 대상 그림의 ACC는 무차원, ubRMSE는 m³/m³. 정규화 오차성분은 이 그림에 표시하지 않음 |

## 발견 사항과 영향

수치 불일치 여부와 오차는 아래 실행 결과에 기록한다. 월별 지점 구성 변화는 계절 효과와 분리해 해석해야 한다. SMAP 선은 DA_single의 model-space scaled 관측에 대한 공통 baseline이며, 두 run에 실제 투입된 각각의 관측 skill 비교가 아니다. 평가 연도를 climatology 학습에 포함하므로 out-of-sample 검증으로 해석할 수 없다. SMAP skill과 DA 개선의 연관성은 개별 동화 update의 인과 효과를 증명하지 않는다.

모델은 공통 3시간 자료 8개가 모두 있을 때 daily 평균을 만들지만, ISMN은 QC good hourly 평균의 가용값으로 daily 평균을 만든다. 따라서 일자는 같아도 하루 내 관측 시각과 coverage까지 동일하다는 보장은 없다. 이 감사는 cache builder 코드를 읽었으며 원 ISMN hourly coverage를 다시 검산하지 않았다.

## 재현 및 보존

저장소 루트에서 다음 두 명령을 순서대로 실행한다.

```bash
/home/gychoi/anaconda3/envs/kite/bin/python docs/audits/smap_skill_20261019/audit_comparison.py
/home/gychoi/anaconda3/envs/kite/bin/python docs/audits/smap_skill_20261019/check_samples.py
```

기존 점검 스크립트를 별도 작업 경로로 복사하여 현재 입력에 다시 실행했다. 대상 함수 재호출 없이 별도 climatology 및 상관 계산으로 daily cache에서 재계산했다. cache 이전 L2 QC, 모델 실행, 실제 동화 배치는 검증 범위 밖이다. notebook 자체를 실행하면 원래 출력에 쓰기 때문에 실행하지 않았다. 대신 notebook과 Python 코드의 일치를 확인하고 계산을 별도로 재현했다. 입력 SHA-256을 실행 전후 비교했다. 수치 허용 오차는 1e-10(ACC 무차원, ubRMSE m³/m³, TE (m³/m³)²); 차이 열 검산은 1e-12이다.

## 이번 실행 결과

두 스크립트 모두 exit 0. notebook과 대응 Python 코드가 일치했다. 356개 지점(Filtered 312개), 관측일 지점·연도·월 20,401개와 최종 집계 96행을 재현했다. 평가 축은 2015-01-02–2025-12-31이며, 공통 anomaly 관측일의 전체 범위는 2015-03-31–2025-12-30이다. 월별 채택 표본은 2015–2025년, 월당 8–31일이다. 지점마다 시작·끝과 가용 연도는 다르다.

| 재계산 검사 | 최대 절대 차이 |
|---|---:|
| 월별 ACC | 8.066e-14 |
| 월별 ubRMSE | 9.975e-17 m³/m³ |
| TE−(VE+CE) | 1.388e-17 (m³/m³)² |
| Screening margin | 9.714e-17 |
| 최종 평균 ACC/ubRMSE | 3.331e-16 |
| 최종 ΔACC | 1.527e-16 |
| 최종 ΔubRMSE | 1.059e-16 m³/m³ |

오차는 모두 허용 범위 이내이다. 월별 N, screening, class 구성도 통과했다. 검사 입력의 실행 전후 SHA-256이 동일했다. 수치 수정은 필요하지 않다.

All은 월별 173–349개, filtered는 153–306개 지점이다. 12개 calendar month 모두에 한 번 이상 기여한 공통 지점은 각각 150/135개이다. 이는 매년 모든 달이 완비된 지점 수를 뜻하지 않는다. 월간 구성은 monthly_site_composition.csv에 key 단위로 저장했다.

권장되는 설명 보완은 (1) 월별 표본 변화, (2) SMAP이 DA_single-scaled 기준이라는 점, (3) 공통 평가 날짜와 climatology 학습 날짜의 차이를 명시하는 것이다. 고정 지점의 계절 변화를 주장하려면 별도 민감도 분석이 필요하지만, 이번 요청에서는 기존 모집단과 notebook을 변경하지 않았다.

산출물: audit.json(입력·hash·재계산 오차), supplement.json(코드 일치·날짜·지점 교집합·차이 검산), recalculated_station_year_month.csv(독립 월별 값), date_support.csv(지점 날짜), sample_counts.csv(월별 N 범위), monthly_site_composition.csv(월별 key와 교체), 두 재현 스크립트.
