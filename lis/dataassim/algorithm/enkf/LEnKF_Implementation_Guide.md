# LEnKF 구현 가이드

## 목표

LISF의 현재 1D-Var 방식을 **진짜 LEnKF**로 변경하여:
- 각 격자점이 **주변 관측 모두** 사용
- **Localization radius** 내 관측만 선택
- **Gaspari-Cohn 함수**로 거리 가중치 적용

---

## 수정 필요 파일

### 1. `lis/core/LIS_lsmMod.F90`
**함수:** `LIS_lsm_DAmapTileSpaceToObsSpace` (라인 1754-1809)

**현재:**
```fortran
st_id = gid  ! 하나의 observation grid
en_id = gid
```

**수정 후:**
```fortran
! localization_radius 내 모든 관측 검색
do i = 1, total_obs
    if (distance(tile, obs(i)) < radius) then
        if (st_id < 0) st_id = i
        en_id = i
    endif
enddo
```

---

### 2. `lis/dataassim/algorithm/enkf/enkf_Mod.F90`
**위치:** 라인 173 (enkf_setup)

**추가:** Localization radius 설정 읽기

```fortran
! 현재 (라인 ~120)
call ESMF_ConfigGetAttribute(LIS_config, enkf_struc(i,k)%localization_factor, &
     label="DA localization radius factor:", rc=status)

! 추가 필요
call ESMF_ConfigGetAttribute(LIS_config, enkf_struc(i,k)%localization_radius_km, &
     label="DA localization radius (km):", rc=status)
if(status.ne.0) then
    enkf_struc(i,k)%localization_radius_km = 200.0  ! default 200km
endif
```

---

### 3. `lis/dataassim/algorithm/enkf/enkf_types.F90`
**추가:** Localization radius 변수

```fortran
type enkf_dec_type
    ! 기존 변수들...
    real :: localization_factor

    ! 추가
    real :: localization_radius_km  ! Localization radius in km
end type enkf_dec_type
```

---

## 단계별 구현

### Phase 1: 준비 (1시간)

**1.1 백업**
```bash
cd /land1/user/gychoi/LIS/test_merge_DA_LISF/lis/core
cp LIS_lsmMod.F90 LIS_lsmMod.F90.bak_before_LEnKF

cd /land1/user/gychoi/LIS/test_merge_DA_LISF/lis/dataassim/algorithm/enkf
cp enkf_Mod.F90 enkf_Mod.F90.bak_before_LEnKF
cp enkf_types.F90 enkf_types.F90.bak_before_LEnKF
```

**1.2 관측 데이터 구조 파악**

중요: Observations 배열이 어떻게 구성되어 있는지 확인 필요
```fortran
! enkf_Mod.F90에서
type(obs_type) :: Observations(Nobs)

! obs_type 구조
type obs_type
    integer :: species
    real    :: lon, lat
    real    :: value
    real    :: std
    logical :: assim
end type
```

---

### Phase 2: 핵심 수정 (2-3시간)

**2.1 LIS_lsmMod.F90 수정**

**위치:** 라인 1754-1809

**수정 방법:**

```fortran
subroutine LIS_lsm_DAmapTileSpaceToObsSpace(n,k,tileid,st_id,en_id)
    use LIS_coreMod
    use LIS_DAobservationsMod

    ! ... 기존 선언 ...

    ! 추가 변수
    real    :: tile_lat, tile_lon
    real    :: localization_radius_km
    integer :: i, obs_count
    real    :: dist_km

    ! 타일 좌표
    tile_lat = LIS_domain(n)%grid(...)%lat
    tile_lon = LIS_domain(n)%grid(...)%lon

    ! Localization radius (임시로 하드코딩, 나중에 config에서 읽기)
    localization_radius_km = 200.0

    ! 관측 검색
    st_id = -1
    en_id = -1
    obs_count = 0

    ! CRITICAL: Observations 배열 접근 방법 확인 필요!
    ! 현재는 의사코드:

    ! 방법 1: LIS_obs_domain을 통해 검색
    do r = 1, LIS_rc%obs_lnr(k)
        do c = 1, LIS_rc%obs_lnc(k)
            gid = LIS_obs_domain(n,k)%gindex(c,r)

            if (gid > 0) then
                ! observation grid cell의 좌표
                call ij_to_latlon(LIS_obs_domain(n,k)%lisproj, &
                     real(c), real(r), obs_lat, obs_lon)

                ! 거리 계산
                dist_km = haversine_distance(tile_lat, tile_lon, &
                                             obs_lat, obs_lon)

                if (dist_km <= localization_radius_km) then
                    if (st_id < 0) st_id = gid
                    en_id = gid
                    obs_count = obs_count + 1
                endif
            endif
        enddo
    enddo

    ! 디버그
    if (mod(tileid, 1000) == 1) then
        write(LIS_logunit,*) '[LEnKF] Tile ', tileid, ' found ', &
             obs_count, ' obs within ', localization_radius_km, ' km'
    endif

end subroutine
```

**2.2 거리 계산 함수 추가**

`LIS_lsmMod.F90` 또는 별도 모듈에 추가:

```fortran
real function haversine_distance(lat1, lon1, lat2, lon2) result(dist_km)
    implicit none
    real, intent(in) :: lat1, lon1, lat2, lon2
    real, parameter :: R_earth = 6371.0  ! km
    real, parameter :: deg2rad = 3.141592653589793 / 180.0
    real :: dlat, dlon, a, c

    dlat = (lat2 - lat1) * deg2rad
    dlon = (lon2 - lon1) * deg2rad

    a = sin(dlat/2.0)**2 + &
        cos(lat1*deg2rad) * cos(lat2*deg2rad) * sin(dlon/2.0)**2
    c = 2.0 * atan2(sqrt(a), sqrt(1.0-a))

    dist_km = R_earth * c
end function haversine_distance
```

---

### Phase 3: 설정 및 테스트 (1-2시간)

**3.1 설정 파일 수정**

`lis.config`에 추가:
```
DA localization radius (km): 200.0
```

**3.2 컴파일**
```bash
cd /land1/user/gychoi/LIS/test_merge_DA_LISF
./compile
```

**3.3 간단한 테스트**
```bash
# 짧은 기간만 실행
# 로그에서 확인:
grep "LEnKF.*found" lislog.0000
```

**예상 출력:**
```
[LEnKF] Tile 1 found 15 obs within 200.0 km
[LEnKF] Tile 1001 found 8 obs within 200.0 km
```

---

## 구현 시 주의사항

### ⚠️ Critical Issues

**1. Observations 배열 구조 확인 필요**

현재 코드에서:
```fortran
! enkf_Mod.F90:318
allocate(Observations(Nobs))

! enkf_Mod.F90:320
call generateObservations(n, k, Nobjs, Nobs, LIS_OBS_State(n,k), &
     LIS_OBS_Pert_State(n,k), Observations)
```

**문제:**
- `Observations` 배열이 어떻게 구성되어 있는지?
- 각 원소의 `lon, lat`이 정확한 관측 위치인지?
- Observation grid cell의 대표 좌표인지?

**해결:**
1. `generateObservations` 함수 확인
2. `Observations(:)%lon, lat` 값 디버그 출력
3. 실제 관측 위치인지 확인

---

**2. 관측 순서 보장**

수정된 함수가 반환하는 `st_id`, `en_id`는:
- **연속된 인덱스**여야 함
- `Observations(st_id:en_id)` 배열 접근 가능해야 함

**현재 문제:**
- Observation grid를 순회하면서 선택하면 인덱스가 불연속적
- `gid`는 grid index이지 Observations 배열 인덱스가 아닐 수 있음

**해결 방안 1:** Observation list 재구성
```fortran
! 임시 배열에 선택된 관측 저장
allocate(selected_obs_idx(max_obs))
count = 0
do i = 1, total_obs
    if (within_radius(i)) then
        count = count + 1
        selected_obs_idx(count) = i
    endif
enddo

! 이후 enkf_Mod.F90에서 이 리스트 사용
```

**해결 방안 2:** Observation 배열을 지역적으로 재구성
- 더 복잡하지만 메모리 효율적

---

**3. 메모리 관리**

LEnKF에서:
```fortran
N_selected_obs = obs_count * Nobjs
```

만약 200km radius 내에 평균 50개 관측:
- `N_selected_obs = 50 * 1 = 50`
- `Obs_cov`: 50 × 50 행렬
- 타일 수: ~수만 개

**메모리:** 관리 가능 (현재 1개 → 50개로 증가)

---

**4. Gaspari-Cohn 함수 확인**

`enkf_general.F90`의 `get_gaspari_cohn` 함수:
```fortran
dweight = get_gaspari_cohn(dx, dy, xcompact, ycompact)
```

이 함수가:
- 제대로 작동하는지
- xcompact, ycompact 단위가 degree인지 km인지 확인

---

## 검증 방법

### 1. Localization 효과 확인

**실험:**
```python
# LEnKF with radius=50km
incr_50km = xr.open_dataset('.../incr.nc')

# LEnKF with radius=200km
incr_200km = xr.open_dataset('.../incr.nc')

# LEnKF with radius=1000km (거의 EnKF)
incr_1000km = xr.open_dataset('.../incr.nc')

diff_50_200 = incr_50km - incr_200km
print(diff_50_200.max())  # 예상: 0.001~0.01 (유의미한 차이)
```

**기대:**
- Radius 변화에 따라 increment 차이
- Radius ↑ → 더 많은 관측 영향 → increment ↑

---

### 2. 공간 패턴 확인

```python
# 관측 위치
obs_lons = [...]
obs_lats = [...]

# Increment 시각화
plt.figure(figsize=(15, 5))

plt.subplot(131)
incr_1d_var.plot(vmin=-0.01, vmax=0.01)
plt.scatter(obs_lons, obs_lats, c='red', s=5)
plt.title('Current (1D-Var)')

plt.subplot(132)
incr_lenkf.plot(vmin=-0.01, vmax=0.01)
plt.scatter(obs_lons, obs_lats, c='red', s=5)
plt.title('LEnKF')

plt.subplot(133)
(incr_lenkf - incr_1d_var).plot()
plt.scatter(obs_lons, obs_lats, c='red', s=5)
plt.title('Difference')
```

**기대:**
- LEnKF: 관측 주변에 부드러운 영향 확산
- 1D-Var: 관측 위치에서만 불연속적 업데이트

---

## 예상 결과

### Before (1D-Var)
```
Tile 1: 1 obs
Tile 2: 1 obs (or 0 obs)
Tile 3: 1 obs
...
Total obs used per tile: 0~1
Spatial pattern: 점점이 떨어진 업데이트
```

### After (LEnKF with 200km radius)
```
Tile 1: 15 obs within 200km
Tile 2: 8 obs within 200km
Tile 3: 23 obs within 200km
...
Total obs used per tile: 0~50
Spatial pattern: 부드러운 공간적 연속성
```

### Increment 차이
```
1D-Var:  Max increment ~ 0.005 at obs points only
LEnKF:   Max increment ~ 0.010 with smooth spatial pattern
Difference: 0.001~0.005 (significant!)
```

---

## 타임라인

| Phase | Task | Time | Output |
|-------|------|------|--------|
| 1 | 백업 및 구조 파악 | 1h | 백업 파일 |
| 2 | 코드 수정 | 2-3h | 수정된 .F90 |
| 3 | 컴파일 및 디버그 | 1-2h | 실행 파일 |
| 4 | 짧은 테스트 | 1h | 로그 확인 |
| 5 | 전체 실행 | 6-12h | 결과 비교 |
| **Total** | | **11-19h** | **LEnKF 구현** |

---

## 다음 단계

**우선순위 1: 관측 구조 파악** ⭐⭐⭐
```bash
# generateObservations 함수 찾기
grep -r "subroutine generateObservations" lis/
```

**우선순위 2: 간단한 프로토타입**
- 먼저 작은 도메인에서 테스트
- 로그로 관측 개수 확인

**우선순위 3: 본격 구현**
- LIS_lsmMod.F90 수정
- 컴파일 및 테스트

---

**작성일:** 2025-12-30
**버전:** LISF 7.5
**목표:** 1D-Var → LEnKF
