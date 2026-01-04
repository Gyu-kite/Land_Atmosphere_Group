# EnKF 모듈 비교 분석: enkf_Mod.F90 vs enkf_Mod.F90.bak

## 요약

현재 버전(`enkf_Mod.F90`)은 EnKF 알고리즘에 **사용자 설정 가능한 공간 지역화(spatial localization)** 기능을 추가하여, 기존의 하드코딩된 지역화 파라미터를 설정 파일에서 조정 가능하도록 개선했습니다.

---

## 상세 변경 사항

### 1. **데이터 구조 확장** (73번째 줄)

**현재 버전 (`enkf_Mod.F90`):**
```fortran
type, public ::  enkf_dec
   logical     :: fileOpen
   real, allocatable :: innov(:)
   real, allocatable :: forecast_var(:) !HPHt
   real, allocatable :: anlys_res(:)
   real, allocatable :: anlys_incr(:,:)
   real, allocatable :: norm_innov(:)
   real, allocatable :: k_gain(:,:)
   real :: localization_factor = 5.0  ! for localization  ← 새로 추가된 필드
end type enkf_dec
```

**백업 버전 (`enkf_Mod.F90.bak`):**
```fortran
type, public ::  enkf_dec
   logical     :: fileOpen
   real, allocatable :: innov(:)
   real, allocatable :: forecast_var(:) !HPHt
   real, allocatable :: anlys_res(:)
   real, allocatable :: anlys_incr(:,:)
   real, allocatable :: norm_innov(:)
   real, allocatable :: k_gain(:,:)
end type enkf_dec
```

**영향:** 지역화 반경 승수를 저장하는 새 필드가 추가되었으며, 기본값은 5.0입니다.

---

### 2. **설정 파일 읽기 기능 추가** (154-167번째 줄)

**현재 버전 (`enkf_Mod.F90`):**
```fortran
!----------------------------------------------------------------------------
! Read localization parameters from config file
!----------------------------------------------------------------------------
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             enkf_struc(n,k)%localization_factor, &
             label="EnKF localization radius factor:", rc=status)
        if(status.ne.0) then
           enkf_struc(n,k)%localization_factor = 5.0  ! 기본값
           write(LIS_logunit,*) '[INFO] EnKF localization radius factor not found in config, using default: 5.0'
        else
           write(LIS_logunit,*) '[INFO] EnKF localization radius factor set to:', &
                enkf_struc(n,k)%localization_factor
        endif
     enddo
```

**백업 버전 (`enkf_Mod.F90.bak`):**
```fortran
! (이 섹션이 존재하지 않음)
```

**영향:**
- LIS 설정 파일에서 `EnKF localization radius factor:` 값을 읽어옴
- 설정 파일에 지정되지 않은 경우 기본값 5.0 사용
- 설정된 값 또는 기본값을 로그에 기록하여 사용자가 확인 가능

---

### 3. **지역화 반경 계산 방식 변경** (275-276번째 줄 vs 293-299번째 줄)

**현재 버전 (`enkf_Mod.F90`):**
```fortran
call LIS_getDomainResolutions(n,dx,dy)
state_size = LIS_surfaceModel_DAgetStateSpaceSize(n,k)

!xcompact = dx*10.0
!ycompact = dy*10.0
xcompact = dx * enkf_struc(n,k)%localization_factor
ycompact = dy * enkf_struc(n,k)%localization_factor

write(LIS_logunit,*) '[INFO] EnKF localization radius (degrees): ', &
    sqrt(xcompact**2 + ycompact**2)
```

**백업 버전 (`enkf_Mod.F90.bak`):**
```fortran
call LIS_getDomainResolutions(n,dx,dy)
state_size = LIS_surfaceModel_DAgetStateSpaceSize(n,k)

xcompact = dx*10.0
ycompact = dy*10.0
```

**영향:**
- **이전:** x, y 방향 모두 10.0으로 하드코딩된 승수 사용
- **현재:** 설정 파일 또는 기본값(5.0)으로부터 `localization_factor` 사용
- **현재:** 실제 지역화 반경(도 단위)을 로그에 출력하여 진단에 활용 가능
- 주석 처리된 라인에 이전의 하드코딩 방식이 보존됨

---

## 설정 파일 사용법

새로운 기능을 사용하려면 LIS 설정 파일에 다음 라인을 추가하세요:

```
EnKF localization radius factor: 5.0
```

**유효한 값:** 양의 실수
- **낮은 값** (예: 3.0): 작은 지역화 반경, 공간적으로 더 독립적
- **높은 값** (예: 10.0): 큰 지역화 반경, 공간적으로 더 많이 평활화
- **기본값:** 5.0 (설정되지 않은 경우)

---

## 기능적 영향

### 공간 지역화 (Spatial Localization)
지역화 인자(factor)는 EnKF 알고리즘에서 관측값이 상태 업데이트에 영향을 미치는 거리를 제어합니다:

```
지역화 반경 = sqrt((dx * factor)^2 + (dy * factor)^2)
```

여기서:
- `dx`, `dy` = 격자 해상도(도 단위)
- `factor` = 설정 가능한 지역화 반경 인자

### 알고리즘 동작
- **factor = 5.0 (새 기본값):** 기존 하드코딩 값(10.0)보다 보수적인 지역화
- **factor = 10.0:** 기존 하드코딩 동작과 동일
- `enkf_analysis()` 호출 시 `xcompact`와 `ycompact`를 사용하여 공간 지역화에 영향

---

## 마이그레이션 주의사항

**하위 호환성:**
- 설정 파일에 `EnKF localization radius factor:`가 지정되지 **않으면**, 코드는 **5.0**을 기본값으로 사용하며, 이는 기존 하드코딩 값 **10.0**과 **다릅니다**
- 기존 동작을 정확히 유지하려면 설정 파일에 `EnKF localization radius factor: 10.0` 추가 필요

**권장사항:**
- 기본값 변경(5.0 vs 10.0)으로 인한 DA 성능 변화를 기존 실험에서 검토
- 새 기본값(5.0)이 고해상도 영역에서 더 나은 지역화를 제공할 수 있음
- 도메인에 최적화하기 위해 다양한 값으로 실험 수행

---

## 코드 품질 개선사항

1. **한글 주석 보존:** `기본값` - 원 개발자의 주석 유지
2. **사용자 친화적 로깅:** 설정된 값과 기본값을 명확한 INFO 메시지로 표시
3. **자기 문서화:** 주석 처리된 기존 코드가 구현의 진화 과정을 보여줌
4. **유연성:** 사용자가 재컴파일 없이 지역화를 조정 가능

---

## 관련 파일

이 변경사항은 다음 파일들의 업데이트가 필요할 수 있습니다:
- LIS 설정 파일 템플릿/문서
- `enkf_general.F90` (`xcompact`와 `ycompact` 파라미터 사용)
- 지역화 인자 파라미터를 설명하는 사용자 문서

---

## 주요 차이점 요약표

| 항목 | 백업 버전 (`.bak`) | 현재 버전 |
|------|-------------------|----------|
| **데이터 구조** | `localization_factor` 필드 없음 | `real :: localization_factor = 5.0` 추가 |
| **설정 읽기** | 없음 | `enkf_setup`에서 설정 파일 읽기 추가 |
| **지역화 계산** | `dx*10.0`, `dy*10.0` 하드코딩 | `dx * localization_factor` 사용 |
| **기본 승수값** | 10.0 (하드코딩) | 5.0 (설정 가능) |
| **진단 로깅** | 없음 | 지역화 반경 로그 출력 추가 |

---

**분석 날짜:** 2025-12-29
**비교한 파일 버전:**
- 현재: `enkf_Mod.F90` (설정 가능한 지역화)
- 백업: `enkf_Mod.F90.bak` (하드코딩된 지역화)
