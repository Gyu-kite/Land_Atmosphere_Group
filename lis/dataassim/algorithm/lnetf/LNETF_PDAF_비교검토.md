# LNETF 구현 vs PDAF 비교 검토

## 1. 알고리즘 단계별 비교

| 단계 | PDAF 구현 | LIS 구현 | 일치 여부 |
|------|-----------|----------|----------|
| 1. Likelihood 계산 | `U_likelihood_l()` 외부 함수 | `compute_likelihood_localized()` 내부 함수 | 방식 상이 |
| 2. 가중치 정규화 | 단순 합 정규화 | **Log-sum-exp trick** | LIS 우수 |
| 3. 변환 행렬 A | `A = diag(w) - w*w^T` | `A = diag(w) - w*w^T` | **동일** |
| 4. 제곱근 계산 | SYEV 고유값 분해 | SSYEV 고유값 분해 | **동일** |
| 5. 변환 적용 | `T = sqrt(N)*sqrt(A)*R + w` | `T = sqrt(N)*sqrt(A)*R + w` | **동일** |
| 6. 앙상블 갱신 | `X^a = X^f * W` | `X^a = X^f * W` | **동일** |

## 2. 핵심 차이점 분석

### 2.1 Likelihood 계산 방식

**PDAF:**
```fortran
! 사용자 정의 함수 호출
CALL U_likelihood_l(domain_p, step, dim_obs_l, obs_l, innov_i, weight)
weights(member) = weight  ! 직접 가중치 반환
```

**LIS:**
```fortran
! Log-likelihood 반환 (수치 안정성)
log_likelihood = -0.5 * sum( w_i * innov_i^2 / R_ii )
```

**평가**: LIS의 log-likelihood 접근법이 수치적으로 더 안정적입니다. 관측이 많거나 innovation이 클 때 underflow를 방지합니다.

### 2.2 Localization 적용 방식

**PDAF:**
```fortran
! Localization은 U_likelihood_l 내부에서 사용자가 구현
! R_localized = R / loc_weight 형태로 적용
```

**LIS:**
```fortran
! 직접 localization 가중치 적용
maha_dist = maha_dist + loc_weights(i) * (innov(i)**2) / Obs_cov(i,i)
```

**수학적 동등성:**
```
PDAF:  L = exp(-0.5 * d^2 / (R/w)) = exp(-0.5 * w * d^2 / R)
LIS:   L = exp(-0.5 * w * d^2 / R)
```
**동등함** - 두 방식 모두 localization 가중치 w가 관측 오차에 적용되어 먼 관측의 영향을 줄입니다.

## 3. 이론적 맹점 분석

### 3.1 관측 오차 공분산 (R) 처리

**현재 구현 (대각선만 사용):**
```fortran
maha_dist = maha_dist + loc_weights(i) * (innov(i)**2) / Obs_cov(i,i)
```

**문제점:**
- R 행렬의 off-diagonal 항 무시
- 상관된 관측 오차 처리 불가

**영향:**
- 대부분의 위성 토양수분 자료는 독립적 관측 오차 가정 → 문제 없음
- 상관된 관측(예: radiance 채널) 사용 시 → 개선 필요

**개선 방안:**
```fortran
! Full Mahalanobis distance 사용
! maha_dist = innov^T * (L * R^{-1} * L) * innov
! where L = diag(loc_weights)
```

### 3.2 Local Analysis Domain 처리

**PDAF:**
- 각 analysis domain에서 local 관측만 사용
- domain_p 인덱스로 병렬화

**LIS:**
- 타일별로 localization radius 내 관측만 필터링
- 동일한 개념으로 구현됨

**평가**: 구현 방식이 다르지만 결과적으로 동일한 local analysis 수행

### 3.3 Weight Degeneracy 처리

**PDAF:**
```fortran
IF (type_winf == 1) THEN
   CALL PDAF_inflate_weights(screen2, dim_ens, limit_winf, weights)
END IF
```

**LIS:**
```fortran
if (type_winf == 1 .and. limit_winf > 0.0) then
   if (n_eff / real(N_ens) < limit_winf) then
      call inflate_weights(N_ens, limit_winf, weights)
   endif
endif
```

**차이점**: LIS는 N_eff/N < limit_winf 조건을 명시적으로 체크, PDAF는 내부에서 체크
**평가**: 동일한 효과

### 3.4 수치 안정성

| 항목 | PDAF | LIS | 비고 |
|------|------|-----|------|
| Likelihood underflow | 직접 exp() | Log-sum-exp trick | LIS 우수 |
| 작은 고유값 처리 | `< 1E-15 → 0` | `< 1E-15 → 0` | 동일 |
| 음수 고유값 | 0으로 설정 | 0으로 설정 | 동일 |
| NaN 체크 | 없음 | 있음 | LIS 우수 |

## 4. 잠재적 문제점 및 권장사항

### 4.1 관측 오차 상관 [낮은 우선순위]

**현재 상태:** 대각 R만 사용
**권장:** 현재 토양수분 DA에는 문제 없음. 향후 다변량 DA 시 개선 필요.

### 4.2 Localization 방식의 이론적 일관성 [확인 완료]

**검토 결과:** LIS 구현이 PDAF의 likelihood-based localization과 수학적으로 동등함.

```
PDAF: w_i = exp(-0.5 * (y-Hx)^T * R^{-1} * (y-Hx))
      with R_localized = R / loc_weight

LIS:  log_w = -0.5 * sum(loc_weight * innov^2 / R)
      → 동등
```

### 4.3 Random Rotation Matrix [확인 필요]

**PDAF:** `rndmat`이 외부에서 전달됨 (global rotation)
**LIS:** 내부에서 생성 (`generate_rndmat`)

**잠재적 문제:**
- PDAF는 모든 local domain에서 동일한 rndmat 사용 (일관성)
- LIS는 각 타일에서 다른 rndmat 생성 가능

**현재 코드 확인:**
```fortran
if (type_trans == 2) then
   ! Deterministic transformation (identity matrix)
   rndmat = 0.0
   do i = 1, N_ens
      rndmat(i,i) = 1.0
   end do
```
**현재 type_trans=2이면 문제 없음** (deterministic)

### 4.4 Ensemble Mean 보존 [확인 완료]

**이론:** W의 각 열의 합 = 1이어야 앙상블 평균 보존

**LIS 구현 확인:**
```fortran
! T = sqrt(N) * sqrt(A) * rndmat + w
! 각 열의 합: sum(T(:,j)) = sqrt(N)*sum(...) + sum(w) = sqrt(N)*0 + 1 = 1
```
- sqrt(A) * rndmat의 각 열 합은 0 (A의 null space 특성)
- sum(w) = 1 (정규화됨)
- 따라서 각 열의 합 = 1 **올바름**

## 5. 종합 평가

### 일치 항목 (정상):
1. 변환 행렬 A = diag(w) - w*w^T
2. 고유값 분해 기반 제곱근 계산
3. 앙상블 변환 X^a = X^f * W
4. Weight inflation 로직
5. Forgetting factor 적용

### LIS 구현이 우수한 항목:
1. **Log-sum-exp trick**: 수치적 안정성 향상
2. **NaN 체크**: 디버깅 용이
3. **상세 진단 출력**: N_eff, loc_weight 등

### 개선 권장 항목:
1. **Off-diagonal R 지원** (향후 다변량 DA 시)
2. **Global random rotation** (type_trans=1 사용 시)

## 6. 결론

**LIS LNETF 구현은 PDAF와 이론적으로 동등하며, 수치 안정성 면에서 일부 개선되어 있습니다.**

현재 출력 결과를 보면:
- 분석 증분이 합리적인 범위 내
- N_eff 진단이 정상 작동
- Localization이 올바르게 적용됨

**정상 작동으로 판단됩니다.**

---
*검토일: 2026-01-02*
*비교 대상: PDAF v2.3 LNETF vs LIS LNETF*
