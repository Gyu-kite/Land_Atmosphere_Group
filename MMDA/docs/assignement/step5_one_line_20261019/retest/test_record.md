# 수정 후 시험 실제 대화 전사

현재 대화에 표시된 사용자 요청과 최종 응답을 수동 전사했다. 원시 세션 로그 내보내기 파일은 아니다. 뒤이어 사용자가 수정 후 Skill 본문도 제공했다. 아래 응답을 기록하면서 보완하거나 재작성하지 않았다.

## 실제 요청문

```text
$audit-da-comparison  README.md를 읽고 아래 notebook의 자료동화 성능 비교가
    동일한 날짜·표본·관측 종류·지표 정의에 기반하는지 점검해 줘.

    notebooks/MMDA_paper_rebuild_20261004/main/new_20261018/
    unified_scaled_daily/monthly_sm_climatology/smap_skill/figures.ipynb
```

## 실제 최종 응답 전체

**각 월 안의 방법 간 비교는 같은 날짜·표본·지표 정의에 기반하며, 수치 오류는 발견되지 않았습니다.** 다만 월별 참여 지점과 climatology 학습 날짜는 동일하지 않습니다.

1. **대상과 비교 기준**

   요청한 `…/monthly_sm_climatology/smap_skill/figures.ipynb`와 상위 계산 코드를 점검했습니다. README는 LNETF·다중센서 구현 개요이므로, 세부 평가 기준은 해당 METHODS 문서와 구현에서 확인했습니다.

   실제 실험은 `/land4/LSM_OUTPUT/MMDA/` 아래 다음 자료입니다.

   - OL: `MP50_OL_init20141231_through202512_direct18c`
   - DA_single: `MP50_EnKF30_nopert_2014_2025`
   - DA_multi: `MP50_39_CLSM_Super_Pooled_anomcdf_2018_2025`

   ISMN을 기준으로 세 실험과 scaled SMAP의 ACC·ubRMSE를 비교합니다. 평가 축은 **2015-01-02–2025-12-31**, 실제 공통 anomaly 관측일 범위는 **2015-03-31–2025-12-30**입니다.

2. **정합성 점검표**

   아래 코드 위치는 `unified_scaled_daily` 기준입니다.

   | 항목 | 판정 | 실제 값·근거 |
   |---|---|---|
   | 방법 간 평가 날짜·표본 | 일치 | `figures.py:530–546`: 다섯 series 공통 유효일, 월당 ≥8일, 비영 anomaly 분산 |
   | 월간 지점 구성 | 불일치 | 월별 All **173–349개**, Filtered **153–306개**. 실제 key 집합도 변화 |
   | 평가 관측 종류 | 일치 | QC-v4 ISMN, SNOTEL 제외. SMAP은 **DA_single DAOBS의 scaled daily mean** |
   | ACC 정의 | 일치 | `figures.py:49–87`: 전 가용 연도 calendar-day 평균을 순환 31일로 평활한 climatology 제거 후 Pearson r |
   | Climatology 학습 날짜 | 불일치 | 모델·ISMN은 네 series 공통일, SMAP은 자체 관측일. 문서에 명시된 차이 |
   | ubRMSE 정의 | 일치 | `figures.py:546`: 같은 평가일의 물리값 오차 `std(ddof=0)` |
   | Screening | 일치 | ≥60일, 비영 분산, 일평균 `r(ISMN,scaled SMAP)−r(ISMN,OL)>−0.2` |
   | 집계·개선 부호 | 일치 | 대상 `figures.py:38–46`: 지점·월 내 연도 평균 → 지점 동일 가중. ΔACC=j−OL, ΔubRMSE=OL−j |
   | 실제 run별 동화 QC·update 표본 | 확인 불가 | 기존 metadata를 사용하며 원본 DAOBS·run 설정 전량은 재검증하지 않음 |

3. **발견 사항과 영향**

   수치 오류와 구분해야 할 **해석 한계**가 있습니다. 월별 곡선 변화에는 계절 효과와 참여 지점 교체가 함께 반영됩니다. SMAP 선은 DA_single-scaled 공통 baseline이므로 두 DA 실험에 각각 투입된 관측의 skill을 비교한 결과가 아닙니다.

   또한 평가 연도가 climatology 학습에 포함됩니다. 모델은 하루 8개 출력을, ISMN은 가용 good 관측을 평균하므로 날짜가 같아도 하루 내 시간 coverage까지 같지는 않습니다.

4. **필요한 수정**

   **수치 수정은 불필요합니다.** 설명에 월별 표본 변화, DA_single-scaled SMAP 출처, climatology 학습일 차이를 명시하는 것이 좋습니다. 원본 notebook과 공유 cache는 변경하지 않았습니다.

5. **재현 정보**

   현재 daily cache에서 **356개 지점(Filtered 312개), 20,401개 지점·연도·월, 최종 96행**을 독립 재계산했습니다. 최대 차이는 ACC **8.07×10⁻¹⁴**, ubRMSE **9.97×10⁻¹⁷ m³/m³**로 허용오차 이내입니다.

   실행 명령은 다음 두 파일을 각각 `kite/bin/python`으로 실행하는 것입니다.

   - `docs/audits/smap_skill_20261019_recheck_01/audit_comparison.py`
   - `docs/audits/smap_skill_20261019_recheck_01/check_samples.py`

   두 실행 모두 성공했습니다. [점검 기록](/land1/user/gychoi/LIS/LISF_MMDA/docs/audits/smap_skill_20261019_recheck_01/TASK_LOG_KO.txt), [수치 결과](/land1/user/gychoi/LIS/LISF_MMDA/docs/audits/smap_skill_20261019_recheck_01/audit.json), [월별 지점 구성](/land1/user/gychoi/LIS/LISF_MMDA/docs/audits/smap_skill_20261019_recheck_01/monthly_site_composition.csv)에 근거를 저장했습니다. 검증 범위는 **daily cache 이후 계산**이며, notebook 자체 재실행과 원자료 QC 검증은 포함하지 않습니다.

## 실행 증거의 범위

기존 코드를 새 경로로 복사하여 이번 입력에서 재실행했다. 최초 재계산 완료 출력 chunk 110528, exit 0; 추가 검산 chunk bff145, exit 0. Skill·reference 읽기 chunk cd7c46. 이 식별자는 대화 출력 식별자이며 파일형 원시 로그가 아니다. 전체 최종 응답은 위 전사이고 commentary 및 모든 도구 호출을 포함한 원시 대화 로그를 추출한 것은 아니다.

실제 실행 명령:

```bash
/home/gychoi/anaconda3/envs/kite/bin/python docs/audits/smap_skill_20261019_recheck_01/audit_comparison.py
/home/gychoi/anaconda3/envs/kite/bin/python docs/audits/smap_skill_20261019_recheck_01/check_samples.py
```

보존 사본은 원래 작업 경로보다 깊으므로 __file__ 기반 저장소 루트 계산을 변경하지 않고 retest/에서 직접 실행할 수 있다고 주장하지 않는다. 위 원래 실행 경로를 사용한다. 이번 기록 단계에서는 시험을 다시 실행하지 않았다.
