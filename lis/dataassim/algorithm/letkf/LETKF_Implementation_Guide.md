# LIS LETKF Implementation Guide

## Local Ensemble Transform Kalman Filter (LETKF) for NASA LISF

**Version**: 1.0
**Date**: January 2026
**Author**: Adapted from PDAF and Hunt et al. (2007)

---

## 1. Overview

### 1.1 LETKF란?

LETKF (Local Ensemble Transform Kalman Filter)는 Hunt et al. (2007)이 제안한 앙상블 칼만 필터의 효율적인 변형입니다. 기존 EnKF와 달리 **앙상블 공간(ensemble space)**에서 분석을 수행하여 계산 효율성을 크게 향상시킵니다.

### 1.2 주요 특징

| 특징 | 설명 |
|------|------|
| **앙상블 공간 분석** | N_ens × N_ens 행렬 연산 (상태 공간보다 작음) |
| **결정론적 필터** | 관측 섭동 불필요 (square-root filter) |
| **자연스러운 병렬화** | 각 격자점에서 독립적으로 분석 가능 |
| **공간 국지화** | Gaspari-Cohn 함수로 원거리 상관 제거 |

### 1.3 다른 필터와 비교

| 항목 | EnKF | LETKF | LNETF |
|------|------|-------|-------|
| 분석 공간 | 상태 공간 | 앙상블 공간 | 앙상블 공간 |
| 관측 섭동 | 필요 | 불필요 | 불필요 |
| 오차 분포 가정 | Gaussian | Gaussian | Non-Gaussian 가능 |
| Weight 계산 | Kalman gain | Kalman gain | Likelihood |
| 계산 복잡도 | O(N_state²) | O(N_ens³) | O(N_ens³) |

---

## 2. LETKF 알고리즘

### 2.1 수학적 정의

Hunt et al. (2007), Physica D, 230, 112-126에 기반합니다.

#### 표기법
- $N$ : 앙상블 크기 (ensemble size)
- $X^f$ : 예보 앙상블 행렬 [N_state × N]
- $\bar{x}^f$ : 예보 앙상블 평균
- $Z = X^f - \bar{x}^f \mathbf{1}^T$ : 예보 앙상블 편차 (perturbation)
- $y^o$ : 관측 벡터
- $H$ : 관측 연산자
- $R$ : 관측 오차 공분산 행렬
- $\rho$ : Forgetting factor (inflation)

#### 핵심 수식

**Step 1: 앙상블 편차 계산**
$$Z = X^f - \bar{x}^f \mathbf{1}^T$$
$$HZ = H(X^f) - \overline{H(X^f)}$$

**Step 2: A⁻¹ 행렬 계산**
$$\tilde{P}^a = A = \left[ \frac{N-1}{\rho} I + (HZ)^T R^{-1} (HZ) \right]^{-1}$$

**Step 3: 가중치 벡터 계산**
$$\bar{w} = A \cdot (HZ)^T R^{-1} (y^o - \overline{H(X^f)})$$

**Step 4: 변환 행렬 계산**
$$W = \sqrt{(N-1) \cdot A} + \bar{w} \mathbf{1}^T$$

**Step 5: 분석 앙상블 생성**
$$X^a = \bar{x}^f \mathbf{1}^T + Z \cdot W$$

### 2.2 알고리즘 흐름도

```
┌─────────────────────────────────────────────────────────────┐
│                    LETKF Analysis                            │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  Input: X^f (forecast ensemble), y^o (observations), R      │
│                                                              │
│  ┌─────────────────────────────────────────────────────┐    │
│  │ Step 1: Compute perturbations                        │    │
│  │   Z = X^f - mean(X^f)                               │    │
│  │   HZ = H(X^f) - mean(H(X^f))                        │    │
│  │   d = y^o - mean(H(X^f))  [innovation]              │    │
│  └─────────────────────────────────────────────────────┘    │
│                          ↓                                   │
│  ┌─────────────────────────────────────────────────────┐    │
│  │ Step 2: Compute A^{-1}                               │    │
│  │   RiHZ = R^{-1} * HZ  [with localization]           │    │
│  │   A^{-1} = (N-1)/ρ * I + HZ^T * RiHZ                │    │
│  └─────────────────────────────────────────────────────┘    │
│                          ↓                                   │
│  ┌─────────────────────────────────────────────────────┐    │
│  │ Step 3: Eigenvalue decomposition                     │    │
│  │   A^{-1} = V * Λ * V^T                              │    │
│  │   A = V * Λ^{-1} * V^T                              │    │
│  └─────────────────────────────────────────────────────┘    │
│                          ↓                                   │
│  ┌─────────────────────────────────────────────────────┐    │
│  │ Step 4: Compute weight vector                        │    │
│  │   w = A * HZ^T * R^{-1} * d                         │    │
│  │     = V * Λ^{-1} * V^T * (HZ^T * R^{-1} * d)        │    │
│  └─────────────────────────────────────────────────────┘    │
│                          ↓                                   │
│  ┌─────────────────────────────────────────────────────┐    │
│  │ Step 5: Compute sqrt(A)                              │    │
│  │   sqrt(A) = V * Λ^{-1/2} * V^T                      │    │
│  │   Scaled: sqrt((N-1)*A) = sqrt(N-1) * sqrt(A)       │    │
│  └─────────────────────────────────────────────────────┘    │
│                          ↓                                   │
│  ┌─────────────────────────────────────────────────────┐    │
│  │ Step 6: Build transform matrix W                     │    │
│  │   W = sqrt((N-1)*A) [* Ω] + w * 1^T                 │    │
│  │   (Ω = optional random rotation)                    │    │
│  └─────────────────────────────────────────────────────┘    │
│                          ↓                                   │
│  ┌─────────────────────────────────────────────────────┐    │
│  │ Step 7: Transform ensemble                           │    │
│  │   X^a = mean(X^f) * 1^T + Z * W                     │    │
│  └─────────────────────────────────────────────────────┘    │
│                                                              │
│  Output: X^a (analysis ensemble)                            │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

---

## 3. Gaspari-Cohn Localization

### 3.1 목적

- **Spurious correlation 제거**: 제한된 앙상블 크기로 인한 원거리 가짜 상관 제거
- **계산 효율성**: 국지적 관측만 사용하여 계산량 감소
- **Rank deficiency 완화**: 효과적인 관측 수 증가

### 3.2 Gaspari-Cohn 함수

5차 다항식 함수로 compact support를 가짐 (Gaspari & Cohn, 1999, Eq. 4.10):

$$C(z) = \begin{cases}
1 - \frac{5}{3}z^2 + \frac{5}{8}z^3 + \frac{1}{2}z^4 - \frac{1}{4}z^5 & 0 \leq z \leq 1 \\
4 - 5z + \frac{5}{3}z^2 + \frac{5}{8}z^3 - \frac{1}{2}z^4 + \frac{1}{12}z^5 - \frac{2}{3z} & 1 < z \leq 2 \\
0 & z > 2
\end{cases}$$

여기서 $z = 2 \cdot |d| / c$, $d$는 거리, $c$는 compact support radius.

### 3.3 LIS 구현 방식

**km 기반 localization (권장):**
```fortran
! 설정: sigma_km = 30 km (Gaspari-Cohn scale)
! compact support = 3 * sigma_km = 90 km
! 검색 반경 = 6 * sigma_km = 180 km

dist_km = haversine_km(tile_lon, tile_lat, obs_lon, obs_lat)
compact_km = 3.0 * sigma_km
d = dist_km / compact_km
loc_weight = gaspari_cohn(d)
```

**격자 기반 localization:**
```fortran
! 설정: factor = 5.0
! xcompact = factor * dx (격자 간격)

dx = tile_lon - obs_lon
dy = tile_lat - obs_lat
d = sqrt((dx/xcompact)^2 + (dy/ycompact)^2)
loc_weight = gaspari_cohn(d)
```

### 3.4 Localization 적용

R⁻¹HZ 계산 시 localization weight를 곱함:

$$[R^{-1}HZ]_{ij} = \frac{\rho_i \cdot [HZ]_{ij}}{R_{ii}}$$

여기서 $\rho_i$는 i번째 관측의 localization weight.

---

## 4. PDAF와의 비교

### 4.1 Forgetting Factor Convention

| 구현 | 수식 | Inflation 의미 |
|------|------|----------------|
| **Hunt et al. (원본)** | $(N-1)/\rho \cdot I$ | ρ < 1 → inflation |
| **PDAF** | forget × (N-1) × I | forget > 1 → inflation |
| **LIS** | $(N-1)/\rho \cdot I$ | ρ < 1 → inflation |

**변환 관계:** `PDAF_forget = 1 / Hunt_rho`

| Inflation 정도 | Hunt/LIS (ρ) | PDAF (forget) |
|---------------|--------------|---------------|
| 없음 | 1.0 | 1.0 |
| 5% | 0.95 | 1.053 |
| 10% | 0.90 | 1.111 |
| 20% | 0.80 | 1.250 |

### 4.2 주요 차이점

| 항목 | PDAF | LIS |
|------|------|-----|
| R 행렬 처리 | Full R 지원 (callback) | 대각 R만 지원 |
| Localization | 사용자 정의 함수 | 내장 Gaspari-Cohn |
| 좌표계 | 추상적 도메인 | 위경도 (Haversine 거리) |
| 후처리 | Column mean 제거 | 없음 (기본 LETKF만) |

### 4.3 코드 비교

**A⁻¹ 계산:**
```fortran
! PDAF (PDAF_letkf_analysis.F90, line 267)
Ainv_l = forget * Ainv_l + tmp_Ainv_l

! LIS (letkf_general.F90, line 261)
Ainv_l = Ainv_l / rho_loc + tmp_Ainv
```

**√A 계산:** (동일)
```fortran
! Both PDAF and LIS
do col = 1, N_ens
   do row = 1, N_ens
      tmp_Ainv(row, col) = V(row, col) / sqrt(eigenvalues(col))
   end do
end do
sqrtA = sqrt(N-1) * tmp_Ainv * V^T
```

---

## 5. LIS 구현 상세

### 5.1 파일 구조

```
lis/dataassim/algorithm/letkf/
├── letkf_Mod.F90        # Main module (6-function interface)
├── letkf_general.F90    # Core analysis routine
├── letkf_types.F90      # Type definitions
└── LETKF_Implementation_Guide.md  # This document
```

### 5.2 모듈 구조

```fortran
module letkf_Mod
  ! 6-function interface for LIS DA framework
  public :: letkf_init        ! 구조체 할당
  public :: letkf_setup       ! 설정 파일 읽기
  public :: letkf_increments  ! 분석 증분 계산 (핵심)
  public :: letkf_update      ! 증분을 상태에 적용
  public :: letkf_diagnostics ! 진단 출력
  public :: letkf_final       ! 메모리 해제
end module

module letkf_general
  ! Core LETKF algorithm
  public :: letkf_analysis    ! 메인 분석 루틴
  public :: gaspari_cohn      ! Localization 함수
  public :: haversine_km      ! 대원 거리 계산
end module

module letkf_types
  ! Data structures
  type :: obs_type           ! 개별 관측 데이터
  type :: obs_param_type     ! 관측 타입별 설정
  type :: update_region_type ! 분석 영역
end module
```

### 5.3 letkf_dec 구조체

```fortran
type :: letkf_dec
   ! Innovation diagnostics
   real, allocatable :: innov(:)         ! y - H(x_mean)
   real, allocatable :: norm_innov(:)    ! normalized innovation
   real, allocatable :: forecast_var(:)  ! HPH^T diagonal
   real, allocatable :: anlys_res(:)     ! analysis residual
   real, allocatable :: anlys_incr(:,:)  ! analysis increment

   ! Localization settings
   real :: localization_factor = 5.0     ! 격자 기반 (기본값)
   real :: localization_scale_km = 30.0  ! km 기반
   logical :: use_km_localization = .false.

   ! LETKF parameters
   real :: forget = 1.0                  ! forgetting factor (ρ)
   integer :: type_trans = 0             ! 0=deterministic, 2=random

   ! Localization diagnostics
   real, allocatable :: n_local_obs(:)      ! 타일별 국지 관측 수
   real, allocatable :: mean_loc_weight(:)  ! 평균 localization weight
   real, allocatable :: max_loc_weight(:)   ! 최대 weight
   real, allocatable :: min_loc_weight(:)   ! 최소 weight
end type
```

---

## 6. 설정 방법

### 6.1 lis.config 옵션

```
# LETKF 알고리즘 선택
Data assimilation algorithm: "LETKF"

# Localization 설정 (둘 중 하나 선택)
# 방법 1: km 기반 (권장, Seo et al. 2021 스타일)
LETKF localization scale (km): 30.0

# 방법 2: 격자 기반
LETKF localization radius factor: 5.0

# Covariance inflation (forgetting factor)
LETKF forgetting factor: 0.95

# Transformation type
# 0 = deterministic (기본값)
# 2 = random rotation
LETKF transformation type: 0
```

### 6.2 권장 설정값

**토양수분 DA (SMAP, ASCAT):**
```
LETKF localization scale (km): 30.0
LETKF forgetting factor: 0.95
LETKF transformation type: 0
```

**참고문헌:**
- Seo, E., et al. (2021): JULES + LETKF with σ=30km
- Tak, Y.-J., et al. (2025): Multi-sensor SMAP+ASCAT LETKF

---

## 7. 출력 진단

### 7.1 Innovation 파일

`LETKF/innov_*.nc`:
- `ninnov_XX`: 정규화된 innovation (|d|/σ)
- `innov_XX`: raw innovation (y - H(x̄))
- `analysis_residual_XX`: 분석 후 잔차
- `forecast_variance_XX`: HPHᵀ 대각 성분

### 7.2 Increment 파일

`LETKF/incr_*.nc`:
- `analysis_increment_XX`: 상태 변수별 분석 증분

### 7.3 로그 출력 예시

```
[INFO] Assimilating Observations using LETKF for DA instance 1
[INFO] LETKF km-based localization:
[INFO]   Sigma (km): 30.0
[INFO]   Compact support (km): 90.0
[INFO] ==========================================
[INFO] LETKF Summary:
[INFO]   Total tiles processed: 10000
[INFO]   Tiles with local observations: 8500
[INFO]   Tiles without observations: 1500
[INFO]   Min local obs per tile: 3
[INFO]   Max local obs per tile: 45
[INFO]   Avg local obs per tile: 12.5
[INFO] ==========================================
```

---

## 8. 수학적 검증

### 8.1 PDAF와의 일치성

| 단계 | PDAF | LIS | 일치 |
|------|------|-----|------|
| Innovation d = y - H(x̄) | ✓ | ✓ | ✅ |
| A⁻¹ 계산 | forget×(N-1) | (N-1)/ρ | ✅ (역수) |
| Eigendecomposition | SYEV | SSYEV | ✅ |
| w = A×HZᵀR⁻¹d | ✓ | ✓ | ✅ |
| √A = √(N-1)×V×Λ^(-½)×Vᵀ | ✓ | ✓ | ✅ |
| W = √A + w | ✓ | ✓ | ✅ |
| Xᵃ = x̄ + Z×W | ✓ | ✓ | ✅ |

### 8.2 대각 R 가정

LIS 구현은 **대각 R 행렬만 지원**합니다:

$$[R^{-1}HZ]_{ij} = \frac{[HZ]_{ij}}{R_{ii}}$$

**적용 가능 조건:**
- 관측 오차가 공간적으로 독립 (대부분의 위성 관측)
- 관측 간 상관이 무시 가능

**제한:**
- 비대각 R (관측 오차 상관)이 필요한 경우 부적합
- 대부분의 토양수분 DA에서는 문제없음

---

## 9. References

1. **Hunt, B.R., Kostelich, E.J. and Szunyogh, I.**, 2007: Efficient data assimilation for spatiotemporal chaos: A local ensemble transform Kalman filter. *Physica D*, 230, 112-126. https://doi.org/10.1016/j.physd.2006.11.008

2. **Gaspari, G. and Cohn, S.E.**, 1999: Construction of correlation functions in two and three dimensions. *Q.J.R. Meteorol. Soc.*, 125, 723-757.

3. **Nerger, L., et al.**, 2005: PDAF - The Parallel Data Assimilation Framework. *Computers & Geosciences*. https://pdaf.awi.de/

4. **Seo, E., et al.**, 2021: Assimilation of SMAP soil moisture retrievals using the Local Ensemble Transform Kalman Filter with JULES. *J. Hydrometeorol.*

5. **Tak, Y.-J., et al.**, 2025: Multi-sensor soil moisture data assimilation using LETKF. (In preparation)

---

## 10. Appendix: Quick Reference

### A. 주요 변수명

| 변수 | 의미 | 차원 |
|------|------|------|
| `N_ens` | 앙상블 크기 | scalar |
| `N_obs` | 관측 수 | scalar |
| `N_state` | 상태 변수 수 | scalar |
| `State_incr` | 상태 앙상블 | [N_state, N_ens] |
| `Obs_pred` | 관측 예측 H(X) | [N_obs, N_ens] |
| `Obs_cov` | 관측 오차 공분산 R | [N_obs, N_obs] |
| `Ainv_l` | A⁻¹ 행렬 | [N_ens, N_ens] |
| `RiHZd` | 가중치 벡터 w | [N_ens] |
| `Asqrt_l` | √A 행렬 | [N_ens, N_ens] |

### B. BLAS/LAPACK 함수

| 함수 | 용도 |
|------|------|
| `SSYEV` | 대칭 행렬 eigendecomposition |
| `SGEMM` | 행렬-행렬 곱 |
| `SGEMV` | 행렬-벡터 곱 |

### C. 디버깅 팁

```fortran
! 첫 번째 격자점에서만 상세 로그 출력
if (gid == 1) then
   write(LIS_logunit,*) '[LETKF] Innovation:', innov_l(1:3)
   write(LIS_logunit,*) '[LETKF] Eigenvalues:', svals(1:3)
   write(LIS_logunit,*) '[LETKF] Weight vector:', RiHZd(1:3)
endif
```
