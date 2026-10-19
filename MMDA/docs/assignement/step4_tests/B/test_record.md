# (b) 대화 전사

현재 대화에서 수동 전사한 실제 요청·최종 답변이다. 원시 세션 로그 내보내기 파일은 아니다. 후속 사용자 메시지로 같은 Skill 본문도 제공되었다.

## 실제 요청

```text
$audit-da-comparison  README.md를 읽고 아래 notebook의 자료동화 성능 비교가
  동일한 날짜·표본·관측 종류·지표 정의에 기반하는지 점검해 줘.

  notebooks/MMDA_paper_rebuild_20261004/main/new_20261018/
  unified_scaled_daily/monthly_sm_climatology/smap_skill/figures.ipynb
```

## 실제 최종 답변

README와 상위 계산 코드·daily cache까지 점검했습니다. **방법 간 비교 날짜·유효 표본·지표 정의는 일치하며, 수치 수정이 필요한 오류는 발견되지 않았습니다.**

| 항목 | 판정 | 확인 결과 |
|---|---|---|
| 날짜·표본 | 일치 | ISMN·OL·두 DA·SMAP의 공통 관측일, 월당 최소 8일 |
| 관측 종류 | 일치 | qc_v4 ISMN과 **DA_single에서 추출한 scaled SMAP** 사용 |
| ACC | 일치 | 모든 가용 연도의 순환 31일 climatology를 제거한 anomaly 상관 |
| ubRMSE | 일치 | 공통 날짜의 일평균 오차 표준편차, `ddof=0` |
| 집계·개선 부호 | 일치 | 지점별 연도 평균 후 지점 동일 가중, 양수가 개선 |
| 월간 지점 구성 | 불일치 | 월마다 기여 지점이 달라짐 |
| 원자료 QC·실제 동화 배치 | 확인 불가 | cache 이전 단계는 재검증하지 않음 |

356개 지점의 **20,401개 지점·연도·월 조합과 최종 96행**을 독립 재계산했습니다. 최대 차이는 ACC `8.07×10⁻¹⁴`, ubRMSE `9.97×10⁻¹⁷`로 허용 오차 이내였습니다.

해석에는 다음 제한이 있습니다.

- 월별 지점 수는 All **173–349개**, filtered **153–306개**입니다. 계절 차이에는 지점 구성 변화가 섞입니다.
- SMAP 선은 공통 DA_single-scaled 기준이며, 각 DA에 실제 투입된 관측을 각각 평가한 결과는 아닙니다.
- 평가 날짜는 같지만 climatology 학습 날짜는 모델/ISMN과 SMAP 사이에 다릅니다. 평가 연도도 학습에 포함됩니다.

원본 notebook·cache는 변경하지 않았습니다. 코드 위치, 재현 명령, 수치 및 월별 지점 구성은 [상세 점검 보고서](/land1/user/gychoi/LIS/LISF_MMDA/docs/audits/smap_skill_20261019/REPORT_KO.md)에 저장했습니다.
