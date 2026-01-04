# LIS CDF 매칭 (Cumulative Distribution Function Matching) 가이드

## 1. 개요

CDF 매칭은 위성 관측 데이터와 모델 시뮬레이션 간의 **통계적 일관성**을 확보하기 위한 bias correction 기법입니다. 관측 데이터의 누적분포함수(CDF)를 모델 데이터의 CDF에 매칭시켜, 두 데이터 간의 체계적인 편차(bias)를 제거합니다.

### 참고 문헌
> Reichle, R. H., and R. D. Koster (2004), **Bias reduction in short records of satellite soil moisture**, *Geophys. Res. Lett.*, 31, L19501, doi:10.1029/2004GL020938.

## 2. 이론적 배경

### 2.1 CDF란?
누적분포함수(CDF)는 확률변수 X가 특정 값 x 이하일 확률을 나타냅니다:

```
F(x) = P(X ≤ x)
```

### 2.2 CDF 매칭 원리

**핵심 아이디어**: 관측값의 CDF 값을 찾고, 동일한 CDF 값에 해당하는 모델 값으로 변환

```
관측값 (obs) → 관측 CDF → CDF 값 (0~1) → 모델 CDF → 변환된 값 (rescaled)
```

**수학적 표현**:
```
obs_rescaled = F_model^(-1)(F_obs(obs_value))
```

여기서:
- `F_obs`: 관측 데이터의 CDF
- `F_model`: 모델 데이터의 CDF
- `F_model^(-1)`: 모델 CDF의 역함수

### 2.3 왜 필요한가?

| 문제점 | 설명 |
|--------|------|
| **Systematic Bias** | 위성 관측과 모델 시뮬레이션 간 체계적 편차 존재 |
| **다른 동적 범위** | 관측과 모델의 값 범위가 다름 (예: 토양수분 0.1~0.3 vs 0.05~0.4) |
| **계절적 편차** | 계절에 따라 편차 패턴이 달라질 수 있음 |
| **지역적 특성** | 지역마다 기후 특성에 따른 편차 차이 |

## 3. LIS에서의 구현

### 3.1 주요 Rescaling 옵션

LIS는 `LIS_dataAssimMod.F90`에서 4가지 rescaling 방법을 제공합니다:

| 옵션 | 함수명 | 설명 |
|------|--------|------|
| **CDF matching** | `LIS_rescale_with_CDF_matching` | 전체 분포를 매칭 (가장 권장) |
| **Normal deviate** | `LIS_rescale_with_normal_deviate_scaling` | 평균/표준편차로 정규화 |
| **Linear scaling** | `LIS_rescale_with_linear_scaling` | 선형 범위 정규화 (0~1) |
| **Anomaly scaling** | `LIS_rescale_with_anomaly` | 아노말리 기반 보정 |

### 3.2 CDF 매칭 핵심 코드

**위치**: `lis/core/LIS_dataAssimMod.F90:636-789`

```fortran
subroutine LIS_rescale_with_CDF_matching(&
     n,             &   ! nest index
     k,             &   ! DA instance index
     nbins,         &   ! CDF bin 개수 (예: 100)
     ntimes,        &   ! 시간 구간 수 (1=전체, 12=월별)
     max_obs_value, &   ! 관측 최대값 (예: 토양수분 0.55)
     min_obs_value, &   ! 관측 최소값 (예: 토양수분 0.02)
     model_xrange,  &   ! 모델 CDF x축 값들
     obs_xrange,    &   ! 관측 CDF x축 값들
     model_cdf,     &   ! 모델 CDF y축 값들 (0~1)
     obs_cdf,       &   ! 관측 CDF y축 값들 (0~1)
     obs_value)         ! 변환할 관측값 (입출력)
```

### 3.3 알고리즘 상세 흐름

```fortran
! Step 1: 현재 월 결정 (월별 CDF 사용 시)
if(ntimes.gt.1) then
   kk = LIS_rc%mo        ! 현재 월
else
   kk = 1                ! 전체 기간 단일 CDF
endif

! Step 2: 각 그리드 포인트에 대해 처리
do t=1,LIS_rc%obs_ngrid(k)

   ! Step 3: IQR (사분위수 범위) 품질 검사
   ! - 관측 IQR이 너무 작으면 CDF가 너무 가파름 → 건너뜀
   index_25 = minloc(abs(obs_cdf(t,kk,:) - 0.25))
   index_75 = minloc(abs(obs_cdf(t,kk,:) - 0.75))
   iqr_obs = obs_xrange(t,kk,index_75(1)) - obs_xrange(t,kk,index_25(1))

   if(iqr_obs .lt. 0.05 * (obs_xrange(t,kk,nbins)-obs_xrange(t,kk,1))) then
      obs_delta(t) = 0   ! 이 그리드는 스킵
   endif

   ! Step 4: 관측값에 해당하는 bin 찾기
   binval = nint((obs_value(col,row)-obs_xrange(t,kk,1)) / obs_delta(t)) + 1
   if(binval.gt.nbins) binval = nbins
   if(binval.le.0) binval = 1

   ! Step 5: 관측 CDF 값 얻기
   cdf_obsval = obs_cdf(t,kk,binval)   ! 0~1 사이 값

   ! Step 6: 모델 CDF에서 동일한 CDF 값에 해당하는 x값 찾기 (역변환)
   i = 1
   do while((model_cdf(t,kk,i) .lt. cdf_obsval) .and. (i.le.nbins))
      i = i + 1
      if(i.gt.nbins) exit
   enddo

   ! Step 7: 변환된 값 할당
   obs_value(col,row) = model_xrange(t,kk,i)

enddo
```

### 3.4 CDF 데이터 구조

```fortran
! 3차원 배열 구조
real :: model_xrange(ngrid, ntimes, nbins)  ! 모델 CDF x축
real :: obs_xrange(ngrid, ntimes, nbins)    ! 관측 CDF x축
real :: model_cdf(ngrid, ntimes, nbins)     ! 모델 CDF y축 (누적 확률)
real :: obs_cdf(ngrid, ntimes, nbins)       ! 관측 CDF y축 (누적 확률)

! 예시 (nbins=10):
! xrange:   [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
! cdf:      [0.05, 0.12, 0.25, 0.40, 0.55, 0.70, 0.82, 0.91, 0.97, 1.00]
```

## 4. CDF 파일 생성 (LDT)

CDF 파일은 **LDT (Land Data Toolkit)**를 사용하여 생성합니다.

### 4.1 LDT 설정 예시

```
# ldt.config
...
DA preprocessing method:              "CDF generation"
DA observation source:                "NASA SMAP soil moisture"
Name of the preprocessed DA file:     ./lis_input.d01.nc

# CDF generation parameters
Number of bins for CDF:               100
Temporal resolution for CDF:          "monthly"  # 또는 "yearly"

# 모델 출력과 관측 데이터 경로
LIS history file:                     ./LIS_HIST/*.nc
Observation data directory:           ./SMAP_L3/
```

### 4.2 CDF 파일 구조 (NetCDF)

```
dimensions:
    ngrid = 15000 ;        # 그리드 포인트 수
    nbins = 100 ;          # bin 개수
    SoilMoist_levels = 1 ; # 토양층 수
    time = 12 ;            # 월별이면 12

variables:
    float SoilMoist_xrange(ngrid, time, SoilMoist_levels, nbins) ;
    float SoilMoist_CDF(ngrid, time, SoilMoist_levels, nbins) ;

attributes:
    :nbins_CDF = 100 ;
    :temporal_resolution_CDF = "monthly" ;
```

## 5. 품질 검사 (Quality Control)

### 5.1 IQR 기반 필터링

LIS는 **사분위수 범위(IQR)**를 사용하여 신뢰할 수 없는 CDF를 필터링합니다:

```fortran
! IQR = 75번째 백분위수 - 25번째 백분위수
iqr_obs = Ub_xrange - Lb_xrange

! IQR이 전체 범위의 5% 미만이면 CDF가 너무 가파름 → 스킵
if(iqr_obs .lt. 0.05 * (obs_xrange(t,kk,nbins)-obs_xrange(t,kk,1))) then
   obs_delta(t) = 0  ! 이 그리드는 CDF 매칭 안 함
endif
```

**이유**: IQR이 작다는 것은:
- 데이터 샘플이 부족하거나
- 데이터가 거의 일정한 값이거나
- CDF가 급격히 증가하여 매칭이 불안정해질 수 있음

### 5.2 경계값 처리

```fortran
! 변환된 값이 허용 범위를 벗어나면 undefined로 설정
if(obs_tmp.gt.max_obs_value) then
   obs_tmp = LIS_rc%udef      ! 최대값 초과 → 제외
endif

if(obs_tmp.le.min_obs_value) then
   obs_tmp = LIS_rc%udef      ! 최소값 미만 → 제외
endif
```

## 6. 사용 예시 (NASA SMAP)

### 6.1 lis.config 설정

```
# Data Assimilation Options
Data assimilation algorithm:                "EnKF"
Number of data assimilation instances:      1

# Observation settings
Data assimilation observation source:       "NASA SMAP soil moisture"
SMAP soil moisture data directory:          ./SMAP_L3/
SMAP soil moisture CDF file:                ./cdf_smapobs.nc

# Rescaling option
Data assimilation observation scaling:      "CDF matching"

# CDF parameters
Number of bins in CDF:                      100
Temporal resolution of CDF:                 "monthly"
```

### 6.2 코드 흐름

```
read_NASASMAPsm.F90:
    ↓
1. SMAP 관측 데이터 읽기
    ↓
2. LIS_readCDFdata() → CDF 파일에서 xrange, cdf 읽기
    ↓
3. LIS_rescale_with_CDF_matching() → 관측값 변환
    ↓
4. 변환된 관측값을 DA 알고리즘에 전달
    ↓
5. EnKF/LNETF 등에서 동화 수행
```

## 7. 다른 Rescaling 방법들

### 7.1 Normal Deviate Scaling

**원리**: 정규분포 가정, 평균과 표준편차로 정규화

```fortran
obs_rescaled = (obs_value - obs_mu) * (model_sigma / obs_sigma) + model_mu
```

**장점**: 계산이 간단, 파라미터 2개만 필요 (평균, 표준편차)
**단점**: 정규분포 가정이 맞지 않으면 부정확

### 7.2 Linear Scaling

**원리**: 최소-최대 범위로 0~1 정규화

```fortran
obs_rescaled = (obs_value - obs_min) / (obs_max - obs_min)
```

**장점**: 가장 간단한 방법
**단점**: 분포 형태를 고려하지 않음

### 7.3 Anomaly Scaling

**원리**: 기후 평균 대비 아노말리만 전달

```fortran
obs_anomaly = obs_value - obs_climatology
obs_rescaled = model_climatology + obs_anomaly
```

**장점**: 기후적 특성 보존
**단점**: 상세한 기후 데이터 필요

## 8. Stratified CDF Matching

### 8.1 개념

**위치 독립적** CDF 매칭으로, 강수량 기후값에 따라 그리드를 그룹화하여 각 그룹별로 CDF를 계산합니다.

```fortran
! 강수 기후값에 따른 bin 결정
strat_binval = nint((target_precip_climo(col,row,kk) - ref_p_min) / delta_p_ref) + 1

! 해당 stratification bin의 CDF 사용
obs_rescaled = model_xrange(strat_binval, kk, matched_bin)
```

### 8.2 장점

- 데이터가 부족한 그리드에서도 안정적인 CDF 추정
- 기후적으로 유사한 지역끼리 그룹화
- CDF Transfer 방식에 활용 (타 센서 데이터 활용)

## 9. 시각화 예시

```
      CDF 매칭 개념도

  1.0 |                    ___○___ Model CDF
      |                ___/   |
      |            ___/       |
  0.7 |--- - - -○-           |  ← 같은 CDF 값
      |        /    \         |
      |       /      \        |
      |      /        \_______|
  0.0 |_____/                 |
      0.1   0.25    0.35   0.5  (값)
            ↑         ↑
         관측값    변환된 값

  동작: 관측값 0.25 → CDF 값 0.7 → 모델에서 0.7에 해당하는 값 0.35로 변환
```

## 10. 문제 해결

### 10.1 일반적인 문제

| 문제 | 원인 | 해결책 |
|------|------|--------|
| 모든 관측값이 undefined | CDF 파일 경로 오류 | 파일 경로 확인 |
| 특정 지역만 undefined | IQR 필터에 걸림 | 데이터 샘플 늘리기 |
| 값이 비현실적으로 큼/작음 | nbins 부족 또는 CDF 품질 | nbins 늘리기, 훈련 기간 연장 |
| 계절적 편차 발생 | ntimes=1 사용 중 | ntimes=12 (월별) 사용 |

### 10.2 디버깅 팁

```fortran
! CDF 매칭 전후 값 출력
write(LIS_logunit,*) 'Before CDF matching: ', obs_value(col,row)
call LIS_rescale_with_CDF_matching(...)
write(LIS_logunit,*) 'After CDF matching: ', obs_value(col,row)
```

## 11. 참고 자료

1. **LIS 소스 코드**: `lis/core/LIS_dataAssimMod.F90`
2. **관측 모듈 예시**: `lis/dataassim/obs/NASA_SMAPsm/read_NASASMAPsm.F90`
3. **LDT CDF 생성**: LDT Users' Guide, Chapter 7
4. **논문**: Reichle and Koster (2004), Bias reduction in short records of satellite soil moisture

---

*작성일: 2026-01-01*
*LIS Version: 7.5*
