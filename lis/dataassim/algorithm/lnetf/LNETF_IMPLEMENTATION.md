# LNETF 구현 상세 문서

**Local Nonlinear Ensemble Transform Filter for LIS**

작성일: 2024년 12월 28일
기반: PDAF (Parallel Data Assimilation Framework)

---

## 목차

1. [개요](#개요)
2. [파일 구조](#파일-구조)
3. [알고리즘 상세](#알고리즘-상세)
4. [각 파일별 상세 설명](#각-파일별-상세-설명)
5. [데이터 흐름](#데이터-흐름)
6. [핵심 수식](#핵심-수식)
7. [사용 방법](#사용-방법)

---

## 개요

### LNETF란?

**LNETF (Local Nonlinear Ensemble Transform Filter)**는 비가우시안 오차 분포와 비선형 관측 연산자를 처리할 수 있는 앙상블 자료동화 방법입니다.

**주요 특징:**
- ✅ **입자 가중치 기반**: 칼만 이득 대신 관측 우도(likelihood)로 가중치 계산
- ✅ **비선형 대응**: 비선형 관측 연산자에 효과적
- ✅ **비가우시안 오차**: 비정규분포 오차에도 적용 가능
- ✅ **Square-root filter**: 수치적 안정성을 위한 제곱근 필터
- ✅ **국지화 지원**: 공간적 국지화(localization) 가능

**참고문헌:**
- Toedter, J., and B. Ahrens, 2015: A Second-Order Exact Ensemble Square Root Filter for Nonlinear Data Assimilation. *Mon. Wea. Rev.*, 143, 1347-1367.
- PDAF documentation: http://pdaf.awi.de/

---

## 파일 구조

### LNETF 디렉토리 구성

```
lnetf/
├── lnetf_Mod.F90          (1,622 lines) - 메인 모듈, LIS 인터페이스
├── lnetf_types.F90        (  104 lines) - 데이터 타입 정의
├── lnetf_general.F90      (  405 lines) - 핵심 LNETF 알고리즘
└── LNETF_IMPLEMENTATION.md            - 본 문서
```

**총 코드 라인: 2,131 lines**

### 파일별 역할

| 파일 | 역할 | 주요 함수 |
|------|------|-----------|
| `lnetf_Mod.F90` | LIS 프레임워크 인터페이스 | `lnetf_init`, `lnetf_setup`, `lnetf_increments`, `lnetf_update` |
| `lnetf_types.F90` | 관측 데이터 구조 정의 | `obs_type`, `obs_param_type` |
| `lnetf_general.F90` | LNETF 핵심 알고리즘 | `lnetf_analysis`, `compute_likelihood_gaussian` |

---

## 알고리즘 상세

### LNETF 알고리즘 흐름도

```
┌─────────────────────────────────────────────────────────────────┐
│ 1. 초기화 (lnetf_init, lnetf_setup)                            │
│    - 구조체 할당                                                 │
│    - 파라미터 읽기 (localization_factor, forget, etc.)          │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ 2. 관측 데이터 준비 (lnetf_increments)                          │
│    - 관측 수집 (Observations)                                   │
│    - 관측 예측 (Obs_pred = H(X^f))                             │
│    - 관측 오차 공분산 (Obs_cov)                                 │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ 3. LNETF 분석 (lnetf_analysis in lnetf_general.F90)            │
│                                                                 │
│    STEP 1: 가중치 계산                                          │
│    ┌──────────────────────────────────────────────────┐        │
│    │ for i = 1 to N_ens:                              │        │
│    │   d_i = y_obs - H(x_i^f)                         │        │
│    │   w_i = exp(-0.5 * d_i^T * R^-1 * d_i)          │        │
│    │ end for                                           │        │
│    └──────────────────────────────────────────────────┘        │
│                              ↓                                  │
│    STEP 2: 가중치 정규화                                        │
│    ┌──────────────────────────────────────────────────┐        │
│    │ w = w / sum(w)                                    │        │
│    │ if sum(w) == 0: w = 1/N_ens (equal weights)      │        │
│    └──────────────────────────────────────────────────┘        │
│                              ↓                                  │
│    STEP 2.5: 유효 샘플 크기 계산                               │
│    ┌──────────────────────────────────────────────────┐        │
│    │ N_eff = 1 / sum(w_i^2)                            │        │
│    └──────────────────────────────────────────────────┘        │
│                              ↓                                  │
│    STEP 3: 변환 행렬 구성                                       │
│    ┌──────────────────────────────────────────────────┐        │
│    │ A = diag(w) - w * w^T                             │        │
│    │   (N_ens × N_ens 행렬)                            │        │
│    └──────────────────────────────────────────────────┘        │
│                              ↓                                  │
│    STEP 4: 고유값 분해                                          │
│    ┌──────────────────────────────────────────────────┐        │
│    │ [V, λ] = eig(A)  (LAPACK dsyev 사용)             │        │
│    │ sqrt(λ_i) 계산                                    │        │
│    │ T = V * diag(sqrt(λ))                             │        │
│    │ sqrt(A) = T * T^T                                 │        │
│    └──────────────────────────────────────────────────┘        │
│                              ↓                                  │
│    STEP 5: 최종 변환 행렬                                       │
│    ┌──────────────────────────────────────────────────┐        │
│    │ W = sqrt(N_ens) * sqrt(A) * R + w                 │        │
│    │   여기서 R은 랜덤 직교 행렬 (현재는 I 사용)       │        │
│    └──────────────────────────────────────────────────┘        │
│                              ↓                                  │
│    STEP 6: 앙상블 변환                                          │
│    ┌──────────────────────────────────────────────────┐        │
│    │ X^a = X^f * W                                     │        │
│    │   (블록 단위로 메모리 효율적 계산)                │        │
│    └──────────────────────────────────────────────────┘        │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ 4. 상태 업데이트 (lnetf_update)                                 │
│    - 분석된 앙상블을 모델 상태에 반영                            │
│    - 진단 정보 출력                                              │
└─────────────────────────────────────────────────────────────────┘
```

---

## 각 파일별 상세 설명

### 1. `lnetf_types.F90` (104 lines)

**목적:** 관측 데이터 구조 정의

#### 주요 타입

##### `obs_type` - 개별 관측

```fortran
type :: obs_type
   integer :: species   ! 관측 타입 식별자 (1=SMAP, 2=ASCAT, etc.)
   integer :: catnum    ! 도메인 내 catchment 번호
   real    :: lon       ! 경도
   real    :: lat       ! 위도
   real    :: value     ! 관측값
   real    :: std       ! 관측 오차 표준편차
   logical :: assim     ! .true. = 동화, .false. = 검증만
end type
```

**사용 위치:**
- `lnetf_general.F90`의 `lnetf_analysis`에서 관측값 추출
- 우도(likelihood) 계산 시 사용

##### `obs_param_type` - 관측 파라미터

```fortran
type :: obs_param_type
   integer :: species         ! 관측 타입
   character(LEN=*) :: path   ! scaling 파라미터 경로
   real :: std                ! 오차 표준편차
   real :: cross_corr         ! 교차 상관계수
   real :: corr_len           ! 상관 길이 [m]
end type
```

**사용 위치:**
- `lnetf_Mod.F90`의 관측 오차 공분산 행렬 구성
- 다중 소스 자료동화 시 species별 처리

---

### 2. `lnetf_Mod.F90` (1,622 lines)

**목적:** LIS 프레임워크와의 인터페이스, 전체 자료동화 워크플로우 관리

#### 주요 데이터 구조

##### `lnetf_dec` 타입

```fortran
type :: lnetf_dec
   logical :: fileOpen

   ! 진단 변수
   real, allocatable :: innov(:)           ! 혁신(innovation)
   real, allocatable :: forecast_var(:)    ! 예측 분산 HPH^T
   real, allocatable :: anlys_res(:)       ! 분석 잔차
   real, allocatable :: anlys_incr(:,:)    ! 분석 증분
   real, allocatable :: norm_innov(:)      ! 정규화된 혁신

   ! LNETF 전용 변수
   real, allocatable :: weights(:)         ! 입자 가중치
   real, allocatable :: transform(:,:)     ! 앙상블 변환 행렬

   ! LNETF 파라미터
   real :: localization_factor = 5.0       ! 국지화 반경 팩터
   integer :: type_forget = 0              ! Forgetting factor 타입
   integer :: type_trans = 0               ! 변환 타입
   integer :: type_winf = 0                ! 가중치 팽창 타입
   real :: forget = 1.0                    ! Forgetting factor
   real :: limit_winf = 0.0                ! 가중치 팽창 한계
   real :: eff_sample_size = 0.0           ! 유효 샘플 크기
end type
```

#### 주요 서브루틴

##### 1) `lnetf_init()` - 초기화

**위치:** Line 101-114
**호출 시점:** LIS 시작 시 1회

**수행 작업:**
```fortran
! 구조체 할당
allocate(lnetf_struc(LIS_rc%nnest, LIS_rc%ndas))
```

##### 2) `lnetf_setup(k)` - 설정

**위치:** Line 121-182
**호출 시점:** 각 DA 인스턴스마다 1회

**수행 작업:**
1. 앙상블 수 검증 (`nensem > 1` 확인)
2. 관측 개수 확인
3. 진단 배열 할당
4. **설정 파일에서 파라미터 읽기:**
   ```fortran
   call ESMF_ConfigGetAttribute(LIS_config, &
        lnetf_struc(n,k)%localization_factor, &
        label="LNETF localization radius factor:", rc=status)
   ```

**설정 파일 예시:**
```
LNETF localization radius factor: 5.0
LNETF forgetting factor: 1.0
LNETF weights inflation limit: 0.8
```

##### 3) `lnetf_increments(n,k)` - 분석 증분 계산

**위치:** Line 190-508
**호출 시점:** 새로운 관측이 있을 때마다

**주요 단계:**

**STEP 1: 관측 데이터 준비**
```fortran
! 관측 수집
call generateObservations(n, k, Nobjs, Nobs, LIS_OBS_State(n,k), &
     LIS_OBS_Pert_State(n,k), Observations)

! 관측 예측 (H(x))
call LIS_surfaceModel_DAGetObsPred(n,k, Obs_pred)

! 관측 섭동
call getObsPert(LIS_OBS_Pert_State(n,k), gsize, N_ens, Nobs, Obs_pert)
```

**STEP 2: 국지화 반경 설정**
```fortran
xcompact = dx * lnetf_struc(n,k)%localization_factor
ycompact = dy * lnetf_struc(n,k)%localization_factor

write(LIS_logunit,*) '[INFO] LNETF localization radius (degrees): ', &
    sqrt(xcompact**2 + ycompact**2)
```

**STEP 3: 각 격자점마다 국지 분석**
```fortran
do i = 1, LIS_rc%ngrid(n)
   ! 국지화 반경 내 관측 선택
   call getSelectedObsNumber(...)

   ! 관측 오차 공분산 행렬 구성
   call assemble_obs_cov(...)

   ! LNETF 분석 수행
   call lnetf_analysis(gid, N_state, N_selected_obs, N_ens, &
        obs_da, obspred_da, obspert_da, Obs_cov, &
        state_incr, state_lon, state_lat, xcompact, ycompact)
end do
```

**STEP 4: 진단 정보 계산**
```fortran
if (LIS_rc%winnov(k).eq.1) then
   ! 혁신, 정규화된 혁신, 예측 분산 계산
   innov = Observations(i)%value - sum(Obs_pred(i,:))/N_ens
   lnetf_struc(n,k)%norm_innov(i) = innov / std_innov
end if
```

##### 4) `lnetf_update(n,k)` - 상태 업데이트

**위치:** Line 516-604
**호출 시점:** `lnetf_increments` 직후

**수행 작업:**
```fortran
! 증분을 모델 상태에 적용
call LIS_surfaceModel_DASetStateVar(n,k)

! 분석 잔차 계산
lnetf_struc(n,k)%anlys_res(i) = Observations(i)%value - &
     sum(Obs_pred(i,:))/N_ens
```

##### 5) `lnetf_diagnostics(n,k)` - 진단 출력

**위치:** Line 613-694
**출력 파일:**
- 혁신 (innovations)
- 정규화된 혁신 (normalized innovations)
- 앙상블 스프레드 (ensemble spread)
- 분석 증분 (analysis increments)

---

### 3. `lnetf_general.F90` (405 lines)

**목적:** LNETF 핵심 알고리즘 구현

#### 주요 서브루틴

##### 1) `lnetf_analysis()` - 핵심 LNETF 분석

**위치:** Line 50-325
**입력:**
- `N_state`: 상태 변수 차원
- `N_obs`: 관측 개수
- `N_ens`: 앙상블 크기
- `Observations`: 관측 데이터 배열
- `Obs_pred`: 관측 예측 H(X^f)
- `Obs_cov`: 관측 오차 공분산 R
- `State_incr`: 예측 앙상블 X^f (출력: 분석 앙상블 X^a)

**출력:**
- `State_incr`: 분석된 앙상블 (변환됨)

**알고리즘 상세:**

**STEP 1: 가중치 계산 (Line 127-149)**
```fortran
! 각 앙상블 멤버에 대해
do n_e = 1, N_ens
   ! 혁신 계산: d = y_obs - H(x_i^f)
   innov_i = obs_values - Obs_pred(:, n_e)

   ! 가우시안 우도 계산
   call compute_likelihood_gaussian(N_obs, innov_i, Obs_cov, weight)

   weights(n_e) = weight
end do
```

**가중치 의미:**
- 관측과 예측이 가까울수록 높은 가중치
- exp(-0.5 * Mahalanobis distance)

**STEP 2: 가중치 정규화 (Line 151-167)**
```fortran
total_weight = sum(weights)

if (total_weight > 1.0e-15) then
   weights = weights / total_weight
else
   ! 모든 가중치가 0인 경우 -> 균등 가중치
   write(LIS_logunit,*) '[WARN] Zero total weight - using equal weights'
   weights = 1.0 / real(N_ens)
endif
```

**STEP 2.5: 유효 샘플 크기 (Line 169-178)**
```fortran
call compute_effective_sample_size(N_ens, weights, n_eff)

! N_eff = 1 / sum(w_i^2)
! N_eff ≈ N_ens: 모든 멤버가 동등하게 기여
! N_eff << N_ens: 소수 멤버만 기여 (가중치 붕괴)
```

**STEP 3: 변환 행렬 A 구성 (Line 183-201)**
```fortran
! A = diag(w) - w * w^T
do j = 1, N_ens
   do i = 1, N_ens
      A_matrix(i,j) = -weights(i) * weights(j)
   end do
end do

do i = 1, N_ens
   A_matrix(i,i) = A_matrix(i,i) + weights(i)
end do
```

**수학적 의미:**
- A는 symmetric positive semi-definite
- rank(A) = N_ens - 1 (하나의 고유값이 0)
- A의 고유벡터는 앙상블 변환의 기저

**STEP 4: 고유값 분해 (Line 203-242)**
```fortran
! LAPACK dsyev 호출
call dsyev('V', 'L', N_ens, A_matrix, N_ens, eigenvalues, work, lwork, info)

! 고유값의 제곱근 계산
do i = 1, N_ens
   if (eigenvalues(i) < 1.0e-15) then
      eigenvalues(i) = 0.0
   else
      eigenvalues(i) = sqrt(eigenvalues(i))
   endif
end do

! T = V * diag(sqrt(λ))
do j = 1, N_ens
   do i = 1, N_ens
      T_matrix(i,j) = A_matrix(i,j) * eigenvalues(j)
   end do
end do

! sqrt(A) = T * T^T
call dgemm('N', 'T', N_ens, N_ens, N_ens, 1.0d0, &
           T_matrix, N_ens, A_matrix, N_ens, 0.0d0, T_tmp, N_ens)
```

**STEP 5: 최종 변환 행렬 (Line 244-272)**
```fortran
! 랜덤 직교 행렬 생성 (현재는 항등 행렬 사용)
rndmat = 0.0
do i = 1, N_ens
   rndmat(i,i) = 1.0
end do

! 스케일링: sqrt(N_ens)
fac = sqrt(real(N_ens))

! W = sqrt(N_ens) * sqrt(A) * R
call dgemm('N', 'N', N_ens, N_ens, N_ens, dble(fac), &
           T_tmp, N_ens, rndmat, N_ens, 0.0d0, T_matrix, N_ens)

! 가중치 추가: W = W + w (브로드캐스트)
do j = 1, N_ens
   do i = 1, N_ens
      T_matrix(i,j) = T_matrix(i,j) + weights(i)
   end do
end do
```

**STEP 6: 앙상블 변환 (Line 278-305)**
```fortran
! 블록 단위 처리로 메모리 절약
maxblksize = 200

do blklower = 1, N_state, maxblksize
   blkupper = min(blklower + maxblksize - 1, N_state)

   ! 예측 앙상블 블록 저장
   do j = 1, N_ens
      ens_blk(1:(blkupper-blklower+1), j) = &
           State_incr(blklower:blkupper, j)
   end do

   ! X^a = X^f * W
   call dgemm('N', 'N', blkupper-blklower+1, N_ens, N_ens, 1.0d0, &
              ens_blk, maxblksize, T_matrix, N_ens, 0.0d0, &
              State_incr(blklower:blkupper,1), N_state)
end do
```

**메모리 효율:**
- 전체 앙상블을 한번에 변환하면 메모리 부족 가능
- 200개 상태변수씩 블록 단위로 처리

##### 2) `compute_likelihood_gaussian()` - 가우시안 우도 계산

**위치:** Line 334-365
**입력:**
- `innov`: 혁신 벡터 (y_obs - H(x))
- `Obs_cov`: 관측 오차 공분산 R

**출력:**
- `likelihood`: exp(-0.5 * d^T * R^{-1} * d)

**구현:**
```fortran
! Mahalanobis 거리 계산 (대각 R 가정)
maha_dist = 0.0
do i = 1, N_obs
   if (Obs_cov(i,i) > 1.0e-15) then
      maha_dist = maha_dist + (innov(i)**2) / Obs_cov(i,i)
   endif
end do

! 가우시안 우도
likelihood = exp(-0.5 * maha_dist)
```

**현재 제한사항:**
- 대각 R만 지원 (관측 오차 독립 가정)
- 향후 개선: Cholesky 분해로 full R 지원

##### 3) `compute_effective_sample_size()` - 유효 샘플 크기

**위치:** Line 374-403
**수식:** N_eff = 1 / Σ(w_i²)

**해석:**
- N_eff = N_ens: 완벽한 균등 가중치
- N_eff = 1: 한 멤버만 유효 (가중치 붕괴)
- N_eff < 0.5*N_ens: 가중치 팽창 고려

---

## 데이터 흐름

### 전체 데이터 흐름도

```
┌─────────────────┐
│ LIS Main Driver │
└────────┬────────┘
         │
         ↓
┌──────────────────────────────────────────────────────┐
│ lnetf_init()                                         │
│ - lnetf_struc 할당                                   │
└────────┬─────────────────────────────────────────────┘
         │
         ↓
┌──────────────────────────────────────────────────────┐
│ lnetf_setup(k)                                       │
│ - 파라미터 읽기 (localization_factor, etc.)          │
│ - 진단 배열 할당                                      │
└────────┬─────────────────────────────────────────────┘
         │
         ↓ (매 동화 사이클)
┌──────────────────────────────────────────────────────┐
│ lnetf_increments(n,k)                                │
│                                                      │
│ ┌────────────────────────────────────────┐          │
│ │ 1. generateObservations()              │          │
│ │    → Observations(1:Nobs)              │          │
│ └────────────────────────────────────────┘          │
│         │                                            │
│         ↓                                            │
│ ┌────────────────────────────────────────┐          │
│ │ 2. DAGetObsPred()                      │          │
│ │    → Obs_pred(Nobs, N_ens)             │          │
│ └────────────────────────────────────────┘          │
│         │                                            │
│         ↓                                            │
│ ┌────────────────────────────────────────┐          │
│ │ 3. assemble_obs_cov()                  │          │
│ │    → Obs_cov(N_obs, N_obs)             │          │
│ └────────────────────────────────────────┘          │
│         │                                            │
│         ↓                                            │
│ ┌────────────────────────────────────────┐          │
│ │ 4. lnetf_analysis() ★                  │          │
│ │                                        │          │
│ │    Input:  X^f(N_state, N_ens)         │          │
│ │            Observations                │          │
│ │            Obs_pred, Obs_cov           │          │
│ │                                        │          │
│ │    Output: X^a(N_state, N_ens)         │          │
│ └────────────────────────────────────────┘          │
│         │                                            │
│         ↓                                            │
│ ┌────────────────────────────────────────┐          │
│ │ 5. 증분 저장                            │          │
│ │    lnetf_struc%anlys_incr = X^a - X^f  │          │
│ └────────────────────────────────────────┘          │
└────────┬─────────────────────────────────────────────┘
         │
         ↓
┌──────────────────────────────────────────────────────┐
│ lnetf_update(n,k)                                    │
│ - 모델 상태에 증분 적용                               │
│ - 분석 잔차 계산                                      │
└────────┬─────────────────────────────────────────────┘
         │
         ↓
┌──────────────────────────────────────────────────────┐
│ lnetf_diagnostics(n,k)                               │
│ - NetCDF 파일 출력                                    │
│   * innovations                                      │
│   * normalized innovations                           │
│   * ensemble spread                                  │
│   * analysis increments                              │
└──────────────────────────────────────────────────────┘
```

### lnetf_analysis 내부 데이터 흐름

```
Input: X^f(N_state, N_ens), Observations(N_obs), Obs_pred(N_obs, N_ens)
  │
  ├─→ STEP 1: 가중치 계산
  │   ┌─────────────────────────────────────┐
  │   │ for i=1:N_ens                       │
  │   │   innov_i = y - H(x_i^f)            │
  │   │   w_i = exp(-0.5*innov^T*R^-1*innov)│
  │   └─────────────────────────────────────┘
  │   weights(N_ens)
  │
  ├─→ STEP 2: 정규화
  │   ┌─────────────────────────────────────┐
  │   │ w = w / sum(w)                      │
  │   └─────────────────────────────────────┘
  │   weights(N_ens) [normalized]
  │
  ├─→ STEP 2.5: N_eff
  │   ┌─────────────────────────────────────┐
  │   │ N_eff = 1 / sum(w^2)                │
  │   └─────────────────────────────────────┘
  │   n_eff (scalar)
  │
  ├─→ STEP 3: 변환 행렬
  │   ┌─────────────────────────────────────┐
  │   │ A = diag(w) - w*w^T                 │
  │   └─────────────────────────────────────┘
  │   A_matrix(N_ens, N_ens)
  │
  ├─→ STEP 4: 고유값 분해
  │   ┌─────────────────────────────────────┐
  │   │ [V, λ] = eig(A)                     │
  │   │ T = V * diag(sqrt(λ))               │
  │   │ sqrt(A) = T * T^T                   │
  │   └─────────────────────────────────────┘
  │   T_tmp(N_ens, N_ens) = sqrt(A)
  │
  ├─→ STEP 5: 최종 변환
  │   ┌─────────────────────────────────────┐
  │   │ W = sqrt(N_ens)*sqrt(A)*R + w       │
  │   └─────────────────────────────────────┘
  │   T_matrix(N_ens, N_ens) = W
  │
  └─→ STEP 6: 앙상블 변환
      ┌─────────────────────────────────────┐
      │ X^a = X^f * W                       │
      │ (블록 단위 BLAS dgemm)               │
      └─────────────────────────────────────┘
      State_incr(N_state, N_ens) ← X^a

Output: X^a(N_state, N_ens)
```

---

## 핵심 수식

### 1. 가중치 계산 (Likelihood)

입자 i에 대한 가우시안 우도:

```
w_i = exp(-0.5 * d_i^T * R^{-1} * d_i)
```

여기서:
- `d_i = y_obs - H(x_i^f)`: 혁신 (innovation)
- `R`: 관측 오차 공분산
- `H(x_i^f)`: 앙상블 멤버 i의 관측 예측

### 2. 가중치 정규화

```
w_i = w_i / Σ(w_j)
```

정규화 후: `Σ(w_i) = 1`

### 3. 유효 샘플 크기

```
N_eff = 1 / Σ(w_i^2)
```

- `N_eff = N_ens`: 균등 가중치
- `N_eff → 1`: 가중치 붕괴

### 4. 변환 행렬

```
A = diag(w) - w * w^T
```

행렬 원소:
```
A_ij = w_i * δ_ij - w_i * w_j
```

여기서 `δ_ij`는 Kronecker delta

### 5. 제곱근 행렬

고유값 분해 `A = V Λ V^T`를 통해:

```
sqrt(A) = V * diag(sqrt(λ_i)) * V^T
```

### 6. 최종 변환

```
W = sqrt(N_ens) * sqrt(A) * R + w
```

여기서:
- `sqrt(N_ens)`: 편향 제거 스케일링
- `R`: 랜덤 직교 행렬 (현재는 I)
- `w`: 가중치 벡터 (각 열에 추가)

### 7. 앙상블 변환

```
X^a = X^f * W
```

행렬 형태:
```
[x_1^a | x_2^a | ... | x_{N_ens}^a] =
[x_1^f | x_2^f | ... | x_{N_ens}^f] * W
```

---

## 사용 방법

### 1. LIS 설정 파일 (lis.config)

```
# 자료동화 알고리즘 선택
Data assimilation algorithm:           LNETF

# LNETF 특화 파라미터
LNETF localization radius factor:      5.0
LNETF forgetting factor:                1.0
LNETF weights inflation type:           0
LNETF weights inflation limit:          0.8
LNETF transformation type:              0

# 일반 DA 설정
Number of ensembles per tile:           20
Assimilation set:                       SMAP soil moisture
```

### 2. 파라미터 설명

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `localization_factor` | 5.0 | 국지화 반경 = dx * factor |
| `forgetting factor` | 1.0 | 공분산 팽창 (1.0 = 팽창 없음) |
| `type_winf` | 0 | 가중치 팽창: 0=없음, 1=N_eff 기반 |
| `limit_winf` | 0.8 | N_eff/N_ens > limit이면 팽창 |
| `type_trans` | 0 | 변환: 0=랜덤 직교, 1=결정론적 |

### 3. 컴파일

```bash
cd /land1/user/gychoi/LIS/test_merge_DA_LISF

# LIS 설정
./configure

# 컴파일
./compile -j 8

# 실행 파일 확인
ls -lh LIS
```

### 4. 실행

```bash
# MPI 실행 (4 프로세스)
mpirun -np 4 ./LIS lis.config

# 출력 확인
tail -f lislog.0000
```

### 5. 출력 파일

LNETF는 다음 진단 파일들을 생성합니다:

```
OUTPUT/LNETF/
├── LNETF.innov.YYYYMMDDHHMM.nc        # 혁신
├── LNETF.spread.YYYYMMDDHHMM.nc       # 앙상블 스프레드
└── LNETF.incr.YYYYMMDDHHMM.nc         # 분석 증분
```

**NetCDF 변수:**
- `ninnov`: 정규화된 혁신
- `innov`: 원시 혁신
- `analysis_residual`: 분석 후 잔차
- `forecast_sigma`: 예측 분산

---

## 향후 개선 사항

### 1. 가중치 팽창 (Weight Inflation)

**현재 상태:** TODO 주석으로 표시
**구현 위치:** `lnetf_general.F90`, Line 180-181

```fortran
! TODO: Implement weights inflation if N_eff/N_ens > limit_winf
if (type_winf == 1) then
   if (n_eff / real(N_ens) < limit_winf) then
      ! 가중치 팽창 수행
      call inflate_weights(weights, n_eff, limit_winf)
   endif
endif
```

### 2. 랜덤 직교 행렬 생성

**현재 상태:** 항등 행렬 사용
**구현 위치:** `lnetf_general.F90`, Line 248-254

```fortran
! TODO: Implement proper random rotation matrix generation
! 현재는 R = I 사용, 향후 PDAF의 generate_rndmat 이식
```

### 3. Full R 행렬 지원

**현재 상태:** 대각 R만 지원
**개선 방법:** Cholesky 분해 사용

```fortran
! Cholesky decomposition: R = L * L^T
call dpotrf('L', N_obs, Obs_cov, N_obs, info)

! Solve: L * y = innov
call dtrsv('L', 'N', 'N', N_obs, Obs_cov, N_obs, innov, 1)

! Mahalanobis distance: ||y||^2
maha_dist = dot_product(innov, innov)
```

### 4. 국지화 구현

**현재 상태:** 준비만 됨
**구현 필요:** Gaspari-Cohn 국지화 함수

```fortran
function get_gaspari_cohn(dx, dy, xcompact, ycompact) result(factor)
   real :: dx, dy, xcompact, ycompact, factor
   real :: dist, ratio

   dist = sqrt(dx**2 + dy**2)
   ratio = dist / sqrt(xcompact**2 + ycompact**2)

   ! Gaspari-Cohn 5차 다항식
   if (ratio < 1.0) then
      factor = 1.0 - 5.0/3.0*ratio**2 + 5.0/8.0*ratio**3 + ...
   else if (ratio < 2.0) then
      factor = 4.0 - 5.0*ratio + 5.0/3.0*ratio**2 + ...
   else
      factor = 0.0
   endif
end function
```

---

## 참고 자료

### 논문

1. **Toedter, J., and B. Ahrens**, 2015: A Second-Order Exact Ensemble Square Root Filter for Nonlinear Data Assimilation. *Mon. Wea. Rev.*, **143**, 1347–1367.
   - LNETF 원본 논문
   - 수학적 유도 및 증명

2. **Tödter, J., and B. Ahrens**, 2015: Generalization of the Ignorance Score: Continuous Ranked Version and Its Decomposition. *Mon. Wea. Rev.*, **143**, 1129–1143.
   - NETF/LNETF 검증 방법

### PDAF 문서

- **PDAF 웹사이트:** http://pdaf.awi.de/
- **PDAF User Guide:** `/home/gychoi/PDAF/doc/PDAF_tutorial.pdf`
- **PDAF Source:** `/home/gychoi/PDAF/src/PDAF_lnetf*.F90`

### LIS 문서

- **LIS Users' Guide:** https://nasa-lis.github.io/LISF/
- **DA Algorithm Guide:** LISF 문서 Chapter 8

---

## 문의 및 기여

**개발자:**
- 기반: PDAF by Lars Nerger & Paul Kirchgessner
- LIS 적용: 2024년 12월 28일

**버그 리포트:**
- LNETF 관련 이슈는 GitHub 또는 LIS support forum

**개선 제안:**
- Pull request 환영

---

**문서 버전:** 1.0
**최종 업데이트:** 2024-12-28
