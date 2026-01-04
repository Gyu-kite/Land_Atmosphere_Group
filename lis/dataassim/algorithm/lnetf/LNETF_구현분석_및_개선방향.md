# LNETF 구현 분석 및 개선 방향

## 1. 개요

이 문서는 현재 LIS(Land Information System)에 구현된 LNETF(Local Nonlinear Ensemble Transform Filter) 알고리즘을 PDAF(Parallel Data Assimilation Framework)의 구현과 비교 분석하여, 현재 구현의 문제점과 개선 방향을 제시합니다.

### 1.1 LNETF 알고리즘 배경

LNETF는 Toedter and Ahrens (2015)에서 제안된 비선형 앙상블 변환 필터입니다:

> J. Toedter and B. Ahrens, "A Second-Order Exact Ensemble Square Root Filter for Nonlinear Data Assimilation", Mon. Wea. Rev. 143 (2015) 1347-1367

**핵심 특징**:
- Kalman gain 대신 **particle weights (likelihood 기반)** 사용
- **비가우시안 오차 분포**에 적합
- **비선형 관측 연산자** 처리 가능
- Transform matrix: `A = diag(w) - w*w^T`
- Eigenvalue decomposition으로 numerical stability 확보

---

## 2. PDAF LNETF 구현 분석

### 2.1 핵심 구조: Local Analysis Domain Loop

PDAF의 LNETF는 **진정한 local filter**로 구현되어 있습니다:

```fortran
! PDAF_lnetf_update.F90:473-640
!$OMP DO schedule(runtime)
localanalysis: DO domain_p = 1, n_domains_p

    ! 1. 현재 local domain의 상태 차원 초기화
    CALL U_init_dim_l(step, domain_p, dim_l)

    ! 2. 전역 앙상블 → 로컬 앙상블 추출
    DO member = 1, dim_ens
        CALL U_g2l_state(step, domain_p, dim_p, ens_p(:, member), dim_l, ens_l(:, member))
    END DO

    ! 3. 로컬 관측 초기화 (핵심!)
    CALL PDAFobs_init_local(domain_p, step, dim_obs_l, ...)

    ! 4. 로컬 분석 수행
    CALL PDAF_lnetf_ana(domain_p, step, dim_l, dim_obs_l, dim_ens, ...)

    ! 5. 로컬 앙상블 → 전역 앙상블 업데이트
    DO member = 1, dim_ens
        CALL U_l2g_state(step, domain_p, dim_l, ens_l(:, member), dim_p, ens_p(:,member))
    END DO

END DO localanalysis
```

### 2.2 Likelihood 계산에서의 Localization

PDAF는 likelihood 계산 시 **거리 기반 가중치**를 적용합니다:

```fortran
! likelihood_l_pdaf.F90 (템플릿)

! Gaspari-Cohn 또는 exponential weight 선택
IF (locweight == 2) THEN
    wtype = 2  ! 5th-order polynomial (Gaspari-Cohn)
END IF

! 각 관측에 대해 거리 기반 가중치 계산
DO i = 1, dim_obs_l
    CALL PDAF_local_weight(wtype, rtype, cradius, sradius, distance_l(i), ...)
END DO

! 가중치 적용된 R^-1 * residual 계산
DO i = 1, dim_obs_l
    Rinvresid_l(i) = ivariance_obs * weight(i) * resid_l(i)
END DO

! 가중된 likelihood 계산
! exp(-0.5 * resid^T * W * R^-1 * resid)
CALL dgemv('t', dim_obs_l, 1, 0.5, resid_l, dim_obs_l, Rinvresid_l, 1, 0.0, likely_l, 1)
likely_l = EXP(-likely_l)
```

### 2.3 PDAF LNETF의 주요 특징

| 구성요소 | PDAF 구현 |
|---------|----------|
| **Local domain loop** | `DO domain_p = 1, n_domains_p` |
| **관측 선택** | `U_init_dim_obs_l` - local domain 주변 관측만 |
| **Localization weight** | `PDAF_local_weight` - Gaspari-Cohn 등 |
| **Likelihood** | 거리 가중 likelihood |
| **Transform matrix** | local domain별 `TA_l` 계산 |
| **앙상블 변환** | `g2l_state` → 분석 → `l2g_state` |

---

## 3. 현재 LIS LNETF 구현의 문제점

### 3.1 핵심 문제: Localization 미구현

`lnetf_general.F90`에서 localization 변수가 선언되지만 **실제로 사용되지 않습니다**:

```fortran
! lnetf_general.F90:70-73 - 선언됨
real, dimension(N_state), intent(in), optional :: State_lon, State_lat
real, intent(in), optional :: xcompact, ycompact

! lnetf_general.F90:127-128 - 플래그 설정됨
apply_localization = (present(State_lon) .and. present(State_lat) &
                     .and. present(xcompact) .and. present(ycompact))

! 하지만 이후 코드에서 apply_localization 변수가 전혀 사용되지 않음!
! 단순히 particle weights 계산 → transform matrix → ensemble 변환
```

### 3.2 EnKF와의 비교

**EnKF (`enkf_general.F90`)의 3D Localization 구현**:

```fortran
! enkf_general.F90:215-253
if (apply_hadamard) then
   ! 3D localization with Hadamard product
   do ii=1,N_state
      do jj=1,N_obs
         ! PHt 계산
         PHt_ij = ...

         ! 거리 계산
         dx = State_lon(ii) - Observations(jj)%lon
         dy = State_lat(ii) - Observations(jj)%lat

         ! Gaspari-Cohn weight 적용 (핵심!)
         dweight = get_gaspari_cohn(dx, dy, xcompact, ycompact)
         PHt_ij = PHt_ij * dweight

         State_incr(ii,n_e) = State_incr(ii,n_e) + PHt_ij*rhs(jj)
      end do
   end do
endif
```

**LNETF (`lnetf_general.F90`)의 현재 구현**:

```fortran
! lnetf_general.F90 - Localization 없음
! 단순히 모든 관측을 동일하게 처리

! Step 1: Likelihood 계산 (localization 없음)
do n_e = 1, N_ens
   innov_i = obs_values - Obs_pred(:, n_e)
   call compute_likelihood_gaussian(N_obs, innov_i, Obs_cov, weight)
   weights(n_e) = weight
end do

! Step 2-6: Transform matrix 계산 및 적용
! ... (localization 고려 없음)
```

### 3.3 문제점 요약

| 항목 | EnKF (3D Local) | LNETF (현재) | PDAF LNETF |
|------|-----------------|--------------|------------|
| Localization 파라미터 | ✅ 있음 | ✅ 있음 | ✅ 있음 |
| 파라미터 전달 | ✅ 전달됨 | ✅ 전달됨 | ✅ 전달됨 |
| 실제 적용 | ✅ Hadamard/GC 적용 | ❌ **미구현** | ✅ 적용됨 |
| 거리 계산 | ✅ `dx, dy` | ❌ 없음 | ✅ `distance_l` |
| 가중치 함수 | ✅ `get_gaspari_cohn` | ❌ 없음 | ✅ `PDAF_local_weight` |
| 분석 단위 | Grid point | Grid point | Local domain |

---

## 4. 개선 방향

### 4.1 방법 1: EnKF 스타일 Localization 추가

EnKF의 Hadamard product 방식을 LNETF likelihood 계산에 적용:

```fortran
! lnetf_general.F90 수정안

subroutine lnetf_analysis(...)
   ...
   ! Likelihood 계산에 localization 적용
   if (apply_localization) then
      do n_e = 1, N_ens
         weighted_likelihood = 0.0

         do i = 1, N_obs
            ! 거리 계산
            dx = State_lon(1) - Observations(i)%lon  ! 상태 위치
            dy = State_lat(1) - Observations(i)%lat

            ! Gaspari-Cohn weight
            loc_weight = get_gaspari_cohn(dx, dy, xcompact, ycompact)

            ! 가중된 innovation
            innov = obs_values(i) - Obs_pred(i, n_e)

            ! 가중된 Mahalanobis distance
            weighted_likelihood = weighted_likelihood + &
                loc_weight * (innov**2) / Obs_cov(i,i)
         end do

         weights(n_e) = exp(-0.5 * weighted_likelihood)
      end do
   else
      ! 기존 비localization 코드
      ...
   endif
   ...
end subroutine
```

### 4.2 방법 2: PDAF 스타일 Local Domain Loop 구현

보다 근본적인 개선으로, PDAF처럼 local analysis domain loop 구현:

```fortran
! lnetf_Mod.F90 수정안

subroutine lnetf_increments(n, k)
   ...
   ! 각 grid point가 아닌 local domain별로 분석
   do domain_p = 1, n_domains_p

      ! 1. Local domain 정의
      call init_local_domain(domain_p, center_lat, center_lon, radius)

      ! 2. 주변 관측 수집 (localization radius 내)
      call collect_local_obs(domain_p, radius, obs_local, n_obs_local)

      ! 3. 로컬 앙상블 추출
      call extract_local_ensemble(domain_p, ens_local)

      ! 4. 거리 기반 likelihood weight 계산
      do i = 1, n_obs_local
         distance = compute_distance(center_lat, center_lon, &
                                     obs_local(i)%lat, obs_local(i)%lon)
         obs_weight(i) = gaspari_cohn(distance, radius)
      end do

      ! 5. 가중된 LNETF 분석
      call lnetf_analysis_local(domain_p, ens_local, obs_local, obs_weight, ...)

      ! 6. 전역 앙상블 업데이트
      call update_global_ensemble(domain_p, ens_local)

   end do
   ...
end subroutine
```

### 4.3 권장 구현 순서

1. **단기 (즉시 적용 가능)**:
   - `lnetf_general.F90`에서 `apply_localization` 플래그 활용
   - EnKF의 `get_gaspari_cohn` 함수 재사용
   - Likelihood 계산에 거리 가중치 적용

2. **중기 (구조적 개선)**:
   - Local observation selection 로직 구현
   - `lnetf_Mod.F90`의 main loop 수정
   - Observation 거리 계산 추가

3. **장기 (PDAF 수준)**:
   - Local domain 개념 도입
   - g2l/l2g state 변환 구현
   - OpenMP 병렬화

---

## 5. 구현 시 고려사항

### 5.1 LNETF vs EnKF Localization의 차이

| 측면 | EnKF | LNETF |
|------|------|-------|
| **적용 대상** | Kalman gain (PHt) | Particle weights (likelihood) |
| **수식** | `PHt * GC_weight` | `exp(-0.5 * d^T * W * R^-1 * d)` |
| **비선형성** | 선형 가정 | 비선형 관측 가능 |

### 5.2 Likelihood에 Localization 적용 방법

PDAF 구현을 참조하면:

```fortran
! Gaussian likelihood with localization
likely_l = 0.0
do i = 1, dim_obs_l
   ! 거리 가중치
   w = gaspari_cohn(distance(i), cradius)

   ! 가중된 R^-1
   Rinv_weighted = w / obs_variance(i)

   ! 가중된 Mahalanobis distance 누적
   likely_l = likely_l + 0.5 * residual(i)**2 * Rinv_weighted
end do

! Likelihood
likely_l = exp(-likely_l)
```

### 5.3 주의사항

1. **Effective Sample Size**: Localization이 강하면 N_eff가 떨어질 수 있음
2. **Weight Inflation**: PDAF는 `type_winf` 옵션으로 weight inflation 지원
3. **Random Rotation Matrix**: 결정론적 vs 확률적 transform 선택 가능

---

## 6. 결론

현재 LIS LNETF 구현은 **Localization이 선언만 되고 실제로 구현되지 않아 사실상 1D 필터**입니다.

EnKF가 3D local filter로 성공적으로 업그레이드된 것처럼, LNETF도 다음이 필요합니다:

1. **Likelihood 계산에 거리 기반 가중치 적용**
2. **Local observation selection 구현**
3. **Gaspari-Cohn 또는 유사 localization function 활용**

PDAF의 `PDAF_lnetf_analysis.F90`과 `likelihood_l_pdaf.F90` 템플릿은 이러한 구현의 좋은 참조가 됩니다.

---

## 참고 파일

| 파일 | 위치 | 설명 |
|------|------|------|
| PDAF LNETF 분석 | `/home/gychoi/PDAF/src/PDAF_lnetf_analysis.F90` | 핵심 분석 루틴 |
| PDAF LNETF 업데이트 | `/home/gychoi/PDAF/src/PDAF_lnetf_update.F90` | Local domain loop |
| PDAF Likelihood 템플릿 | `/home/gychoi/PDAF/templates/.../likelihood_l_pdaf.F90` | Localization weight |
| LIS LNETF 모듈 | `lnetf_Mod.F90` | 현재 구현 |
| LIS LNETF 분석 | `lnetf_general.F90` | 분석 루틴 (localization 미구현) |
| LIS EnKF 분석 | `../enkf/enkf_general.F90` | 3D localization 참조 구현 |

---

## 7. 구현 완료 (2026-01-01)

### 7.1 변경된 파일

**`lnetf_general.F90`** - 핵심 분석 루틴 전면 개선:

1. **Gaspari-Cohn localization 함수 추가**:
   ```fortran
   function gaspari_cohn(d)        ! 5차 다항식 localization
   function get_gaspari_cohn(dx, dy, xcompact, ycompact)  ! 이방성 지원
   ```

2. **새로운 localized likelihood 계산**:
   ```fortran
   subroutine compute_likelihood_localized(N_obs, innov, Obs_cov, loc_weights, likelihood)
   ! likelihood = exp(-0.5 * sum_i( w_i * innov_i^2 / R_ii ))
   ```

3. **Localization 적용 로직**:
   ```fortran
   if (apply_localization) then
      do i = 1, N_obs
         dx = State_lon(1) - Observations(i)%lon
         dy = State_lat(1) - Observations(i)%lat
         loc_weights(i) = get_gaspari_cohn(dx, dy, xcompact, ycompact)
      end do
   endif
   ```

4. **Effective sample size 진단**:
   ```fortran
   subroutine compute_effective_sample_size(N_ens, weights, n_eff)
   ! N_eff = 1 / sum(w_i^2)
   ```

### 7.2 EnKF와의 일관성

| 구성요소 | EnKF | LNETF (신규) |
|---------|------|-------------|
| `gaspari_cohn()` | ✅ | ✅ 동일 |
| `get_gaspari_cohn()` | ✅ | ✅ 동일 |
| 거리 계산 | `dx = State_lon - Obs%lon` | `dx = State_lon(1) - Observations(i)%lon` |
| 가중치 적용 | `PHt * dweight` (Kalman gain) | `loc_weights(i) * innov^2 / R` (likelihood) |

### 7.3 알고리즘 흐름

```
1. Localization 가중치 계산 (if enabled)
   └─ 각 관측에 대해 Gaspari-Cohn 거리 가중치 계산

2. Particle weights (Likelihood) 계산
   └─ compute_likelihood_localized() 사용
   └─ likelihood = exp(-0.5 * sum(loc_weight * innov^2 / R))

3. 가중치 정규화
   └─ weights = weights / sum(weights)

4. Effective sample size 계산
   └─ N_eff = 1 / sum(w^2)

5. Transform matrix 계산
   └─ A = diag(w) - w*w^T
   └─ sqrt(A) via eigenvalue decomposition

6. 앙상블 변환
   └─ X^a = X^f * W
   └─ W = sqrt(N_ens) * sqrt(A) * rndmat + w
```

### 7.4 로그 출력 예시

```
[LNETF] Starting LNETF analysis for grid point 1
[LNETF] N_state= 4  N_obs= 5  N_ens= 12
[LNETF] Localization ENABLED with xcompact= 0.50  ycompact= 0.50
[LNETF] Localization weights (first 3): 1.00  0.85  0.42
[LNETF] Effective sample size: N_eff = 8.5 ( 70.8%)
[LNETF] Analysis complete
```

---

*문서 작성: 2026-01-01*
*구현 완료: 2026-01-01*
*PDAF 버전: 2.3+ (2024-2025)*
*LIS 버전: 7.5*
