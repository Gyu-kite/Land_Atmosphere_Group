# 연구 배경·용어·비교 기준

## 배경과 출처

본 연구는 NASA LISF에서 위성 토양수분을 지면모델에 동화하고 ISMN으로
검증한다. 현재 주 비교는 단일모델 EnKF-NoPert와 pooled MMDA-NoPert이다.
MMDA는 Noah-MP 5.0, Noah 3.9, CLSM의 ensemble을 reference space에서 결합한다.
이 설명은 현재 연구 기록에 근거하며 모든 자료동화 실험에 적용되는 가정은 아니다.

저장소 루트 기준 참고 문서:

- `README.md`: multi-source DA 및 LNETF 중심의 저장소 소개. MMDA의 최신 평가
  정의를 모두 담고 있지는 않으며 2주차 제출 README와 동일한지도 확인되지 않았다.
- `notebooks/MMDA_paper_rebuild_20261004/main/new_20261018/unified_scaled_daily/METHODS_CAPTIONS_KO.md`:
  현재 연구의 공통 날짜, 지표, scaled SMAP 및 screening 정의.
- 위 디렉토리의 `monthly_sm_climatology/smap_skill/METHODS_CAPTIONS_KO.md`:
  월별 공통 관측일 집계 및 baseline 정의.
- `docs/assignement/step2_audit_20261019/RESULTS_KO.md`:
  Skill 추출에 사용한 실제 점검과 확인 범위.

사용자의 현재 요청이 이 기준과 다르면 현재 요청을 우선하고 차이를 밝힌다.
문서와 코드가 다르면 유리한 쪽을 선택하지 말고 불일치로 보고한다.

## 용어와 실제 자료

- OL: 자료동화 없는 Noah-MP 결과.
- DA_single: Noah-MP 30-member EnKF-NoPert. NoPert는 state perturbation OFF이며
  forcing과 observation perturbation까지 없다는 뜻은 아니다.
- DA_multi: 3개 모델 각 10 members의 pooled MMDA-NoPert.
  SURFACEMODEL은 reference MP 평균이며 항상 30-member pooled mean인 것은 아니다.
- ISMN: 현재 연구는 `/land1/ISMN/qc_v4` 기반이며 SNOTEL은 제외한다.
- Scaled SMAP: model-space로 변환된 관측. Native retrieval과 혼동하지 않는다.
  실제 출처가 DA_single인지 DA_multi인지 확인한다. DAOBS 파일 존재만으로
  해당 지점과 시각에 동화되었다고 간주하지 않는다.

## 현재 연구의 지표 기준

공통 daily 날짜를 먼저 맞춘다. 일부 자료가 있는 날짜의 평균과 공통 날짜의
평균을 혼용하지 않는다. UTC daily와 LST diurnal 분석도 구분한다.

ACC는 각 series에서 모든 가용 연도의 calendar-day 평균을 만들고 circular
31-day window로 평활한 climatology를 뺀 anomaly의 Pearson correlation이다.
단순 time-axis rolling이나 leave-one-year-out climatology와 다른 정의이다.
평가 연도가 climatology에 포함되므로 out-of-sample 검증으로 표현하지 않는다.

ubRMSE는 공통 날짜의 absolute daily error에서 평균 error를 제거한 표준편차다.
ddof는 실제 구현에서 확인한다. Anomaly error의 cRMSE와 구분한다.

Centered anomaly error decomposition을 사용할 때:

$$
\mathrm{VE}=(\sigma_x-\sigma_o)^2,\qquad
\mathrm{CE}=2\sigma_x\sigma_o(1-r),\qquad
\mathrm{TE}=\mathrm{VE}+\mathrm{CE}.
$$

동일 날짜·중심화·분산 convention에서만 검산한다. Normalized 지표는
분모가 지점 전체기간 ISMN anomaly variance인지 월별 variance인지 확인한다.

Filtered 기본값은 같은 날짜의 absolute daily Pearson correlation로
`r(ISMN, scaled SMAP) - r(ISMN, OL) > -0.2`이며 최소 60 paired days와
비영 분산을 요구한다. 이 r은 anomaly ACC가 아니다.
All은 상대 skill screening을 적용하지 않은 적격 표본이지 모든 원시 ISMN
지점을 뜻하지 않는다. 동화 지원·기간·공간 조건은 코드에서 별도로 확인한다.

## 월별 그림에서 반복된 점검 포인트

최근 observation-day 그림은 ISMN, OL, DA_single, DA_multi, scaled SMAP이
동일 날짜에 유효하며 월당 최소 8일인 지점·연도·월을 사용한다.
지점/달마다 연도 평균 후 지점을 동일 가중 평균한다.
All-days 그림은 다른 최소 일수 조건을 사용할 수 있으므로 코드를 확인한다.

지점 수가 같아도 구성은 다를 수 있다. 월별 key 집합과 N을 모두 확인하고,
표본 교체와 계절 효과를 구분한다. 고정 지점 sensitivity 분석은 필요성을
제안할 수 있으나 요청 없이 비교 모집단을 바꾸지 않는다.

Delta 지표의 baseline과 부호는 각 그림에서 확인한다. 예를 들어
`ACC_j - ACC_OL`, `ubRMSE_OL - ubRMSE_j`는 모두 양수가 개선이지만,
`TE_multi - TE_baseline`은 음수가 개선이다.

## 검증의 한계

같은 daily cache에서 수치를 재계산하면 cache 이후 계산의 일관성을 검사한 것이다.
L2 reader, QC, domain mapping, 개별 동화 배치까지 확인한 것은 아니다.
정상 판정에는 실제 확인한 코드·파일·날짜 범위를 반드시 붙인다.
