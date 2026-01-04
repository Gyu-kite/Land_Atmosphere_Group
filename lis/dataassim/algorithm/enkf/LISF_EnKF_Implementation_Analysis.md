# LISF EnKF 구현 문제 분석

## 요약

LISF의 현재 EnKF 구현은 **진짜 Ensemble Kalman Filter가 아닙니다**. 각 격자점이 자기 격자의 관측만 사용하는 **1D-Var** 방식에 가깝습니다. Localization 코드는 존재하지만 실질적으로 작동하지 않습니다.

---

## 1. 이론적 EnKF vs LEnKF vs 현재 구현

### 1.1 표준 EnKF (Ensemble Kalman Filter)

**알고리즘:**
```
For each state grid point i:
    st_id = 1
    en_id = total_obs

    obs_selected = ALL observations

    Kalman Gain K = P^f H^T (HPH^T + R)^-1

    x^a_i = x^f_i + K * (y - Hx^f)
```

**특징:**
- ✅ **모든 관측 사용**: 전체 도메인의 모든 관측이 모든 격자점에 영향
- ✅ **공간 정보 활용**: 관측이 먼 곳의 격자점도 업데이트
- ❌ **허위 상관관계**: 앙상블 크기 < 상태변수 개수일 때 문제
- ❌ **계산 비용**: O(N_state × N_obs²)

**수식:**
```
x^a = x^f + P^f H^T (HPH^T + R)^-1 (y - Hx^f)

여기서:
- P^f H^T: 모든 관측과 모든 상태변수 간 공분산
- N_obs × N_state 행렬
```

---

### 1.2 LEnKF (Localized Ensemble Kalman Filter)

**알고리즘:**
```
For each state grid point i:
    st_id = find_first_obs_within_radius(i, radius)
    en_id = find_last_obs_within_radius(i, radius)

    obs_selected = observations within localization radius

    Kalman Gain K = (ρ ⊙ P^f H^T) (ρ ⊙ HPH^T + R)^-1

    x^a_i = x^f_i + K * (y - Hx^f)
```

**특징:**
- ✅ **주변 관측만 사용**: Localization radius 내 관측 선택
- ✅ **Hadamard product (⊙)**: 거리에 따른 가중치 적용
- ✅ **허위 상관관계 제거**: 먼 거리 영향 차단
- ✅ **계산 효율성**: O(N_state × N_local_obs²), N_local_obs << N_obs

**수식:**
```
K_ij = ρ(dist(i,j)) * (P^f H^T)_ij

여기서:
ρ(dist) = Gaspari-Cohn 함수:
- dist < radius: 0 < ρ < 1 (거리에 따라 감소)
- dist > radius: ρ = 0 (영향 없음)
```

---

### 1.3 LISF 현재 구현

**알고리즘:**
```
For each state grid point i:
    gid = map_tile_to_obs_grid(i)

    st_id = gid
    en_id = gid

    obs_selected = observation at grid cell gid ONLY

    x^a_i = x^f_i + K * (y_gid - Hx^f_i)
```

**특징:**
- ❌ **하나의 관측만 사용**: 자기 격자 관측만
- ❌ **공간 정보 없음**: 주변 격자와 독립적
- ❌ **Localization 무의미**: 거리가 항상 ~0
- ✅ **계산 효율성**: O(N_state × 1²) = O(N_state)

**본질:**
- **1D-Var (Optimal Interpolation)에 가까움**
- 각 격자점에서 독립적인 스칼라 업데이트

---

## 2. 코드 분석

### 2.1 문제의 핵심 함수

**파일:** `lis/core/LIS_lsmMod.F90`

**함수:** `LIS_lsm_DAmapTileSpaceToObsSpace`

**위치:** 라인 1754-1809

```fortran
subroutine LIS_lsm_DAmapTileSpaceToObsSpace(n,k,tileid,st_id,en_id)
    ! ...
    ! LSM 타일 → observation grid 매핑

    lat = LIS_domain(n)%grid(...)%lat
    lon = LIS_domain(n)%grid(...)%lon

    ! observation grid에서 대응 cell 찾기
    call latlon_to_ij(LIS_obs_domain(n,k)%lisproj, lat, lon, col, row)
    c = nint(col)
    r = nint(row)

    gid = LIS_obs_domain(n,k)%gindex(c,r)

    ! 문제: 시작과 끝이 같음!
    st_id = gid
    en_id = gid

end subroutine
```

**문제점:**
```fortran
st_id = gid  ! 하나의 관측만
en_id = gid  ! 시작 = 끝
```

→ 하나의 observation grid cell만 반환

---

### 2.2 EnKF 메인 루프

**파일:** `lis/dataassim/algorithm/enkf/enkf_Mod.F90`

**위치:** 라인 365-446

```fortran
do i=1, state_size/LIS_rc%nensem(n)  ! 각 LSM 타일마다

    tileid = (i-1)*LIS_rc%nensem(n)+1

    ! 관측 선택
    call LIS_surfaceModel_DAmapTileSpaceToObsSpace(&
         n, k, tileid, st_id, en_id)

    ! st_id = en_id = gid (하나의 observation grid cell)

    N_selected_obs = (en_id-st_id+1)*Nobjs
    ! = 1 * Nobjs (한 cell의 observation types만)

    ! 관측 추출
    allocate(obs_da(N_selected_obs))
    do kk=st_id,en_id  ! 실제로는 한 번만 실행 (st_id=en_id)
        obs_da(kk) = Observations(kk)
    enddo

    ! EnKF 분석
    call enkf_analysis(gid, N_state, N_selected_obs, N_ens, &
         obs_da, obspred_da, obspert_da, Obs_cov, &
         state_incr, &
         state_lon, state_lat, xcompact, ycompact)  ! Localization 인자

enddo
```

**문제점:**

1. **N_selected_obs = 1 × Nobjs**
   - 하나의 observation grid cell만
   - 주변 관측 미사용

2. **Localization 인자는 전달되지만:**
   - `state_lon, state_lat`: 현재 타일의 좌표
   - `Observations(gid)%lon, lat`: 같은 grid cell의 좌표
   - 거리 ≈ 0 → 가중치 ≈ 1
   - Hadamard product가 의미 없음

---

### 2.3 Localization 코드 (사용되지 않음)

**파일:** `lis/dataassim/algorithm/enkf/enkf_general.F90`

**위치:** 라인 217-249

```fortran
if (apply_hadamard) then
    ! update *with* Hadamard product

    do ii=1,N_state  ! 각 상태변수
        do jj=1,N_obs  ! 각 관측

            ! 거리 계산
            dx = State_lon(ii) - Observations(jj)%lon
            dy = State_lat(ii) - Observations(jj)%lat

            ! Gaspari-Cohn 가중치
            dweight = get_gaspari_cohn(dx, dy, xcompact, ycompact)

            ! 가중치 적용
            PHt_ij = PHt_ij * dweight

        enddo
    enddo
endif
```

**문제점:**

1. **N_obs = 1 (실질적으로)**
   - Loop는 1번만 실행

2. **dx, dy ≈ 0**
   - 같은 grid cell 내
   - dweight ≈ 1

3. **Hadamard product 무의미**
   - 곱하기 1 = 변화 없음

---

## 3. 실험 결과 분석

### 3.1 Localization 유무 비교

**테스트:**
```python
incr_with_loc = xr.open_dataset('SMAPDA_OUTPUT/EnKF/.../incr.nc')
incr_no_loc = xr.open_dataset('SMAPDA_ENKF_OUTPUT/EnKF/.../incr.nc')

diff = incr_with_loc - incr_no_loc
print(diff.max())  # 1.41e-07
print(diff.min())  # -2.08e-07
```

**결과:**
- 차이: **10⁻⁷ m³/m³** (거의 0)
- 정상적인 차이: **10⁻⁴ ~ 10⁻³ m³/m³** 예상

**해석:**
- Localization이 **실질적으로 작동하지 않음**
- 두 실행이 **거의 동일한 결과** 생성

---

### 3.2 디버그 로그 분석

**로그:**
```
[DEBUG-ENKF] Analysis Increment Statistics:
[DEBUG-ENKF]   Analysis increment MAX:  3.3367670E-03
[DEBUG-ENKF]   Analysis increment MIN: -2.5835207E-03
```

**해석:**
- Increment 자체는 정상 (10⁻³ 수준)
- DA는 작동하지만 각 격자점에서 **독립적**으로만 작동

---

## 4. 문제 요약

| 항목 | 표준 EnKF | LEnKF | 현재 LISF |
|------|-----------|-------|-----------|
| **관측 선택** | 모든 관측 | Radius 내 관측 | 자기 격자만 |
| **N_obs per tile** | ~수천 개 | ~수백 개 | **1개** |
| **공간 상관** | 전역 | 국소 | **없음** |
| **Localization** | 불필요 | Hadamard 적용 | **무의미** |
| **계산 복잡도** | O(N×M²) | O(N×m²) | **O(N)** |
| **DA 품질** | 최적 | 준최적 | **제한적** |

**본질:**
- **현재 = 1D-Var (Optimal Interpolation)**
- **NOT = EnKF**
- **NOT = LEnKF**

---

## 5. 수정 방안

### 5.1 진짜 EnKF 구현

**필요한 수정:**

1. **`LIS_lsm_DAmapTileSpaceToObsSpace` 수정**
   ```fortran
   ! 현재
   st_id = gid
   en_id = gid

   ! 수정 후
   st_id = 1               ! 첫 번째 관측
   en_id = total_obs       ! 마지막 관측
   ```

2. **메모리 요구량**
   - 전지구: N_obs ~ 수천 개
   - Obs_cov: N_obs × N_obs 행렬
   - 메모리: **수 GB**

3. **계산 시간**
   - 각 타일마다 모든 관측 처리
   - 시간: **100~1000배 증가**

---

### 5.2 진짜 LEnKF 구현

**필요한 수정:**

1. **`LIS_lsm_DAmapTileSpaceToObsSpace` 수정**
   ```fortran
   ! radius 내 관측 검색
   st_id = -1
   en_id = -1
   count = 0

   do i=1,total_obs
       dist = distance(tile_lat, tile_lon, obs(i)%lat, obs(i)%lon)

       if (dist < localization_radius) then
           if (st_id < 0) st_id = i
           en_id = i
           count = count + 1
       endif
   enddo
   ```

2. **설정 파일에 추가**
   ```
   DA localization radius: 200  # km
   ```

3. **계산 비용**
   - 중간 수준 (EnKF < LEnKF < 현재)
   - 전지구에서 실용적

---

### 5.3 현재 구현 유지 (1D-Var)

**장점:**
- ✅ 매우 빠름
- ✅ 메모리 효율적
- ✅ 병렬화 용이
- ✅ 전지구 규모에서 실용적

**단점:**
- ❌ 공간 정보 미활용
- ❌ 관측 공백 지역 업데이트 없음
- ❌ EnKF/LEnKF라고 부르기 어려움

**권장:**
- 명칭 변경: "Grid-point DA" 또는 "1D-EnKF"
- 문서에 명확히 설명

---

## 6. 결론

### 6.1 현재 상태

**LISF의 "EnKF"는:**
- ❌ 표준 EnKF가 아님
- ❌ LEnKF도 아님
- ✅ 격자점별 독립 DA (1D-Var)

**Localization 코드는:**
- 존재하지만 작동하지 않음
- 각 타일이 하나의 관측만 사용하므로 가중치 계산 무의미

---

### 6.2 권장사항

**옵션 1: 진짜 LEnKF 구현**
- 주변 관측 검색 기능 추가
- 계산 비용 증가하지만 실용 가능
- DA 품질 크게 향상

**옵션 2: 현재 구조 유지**
- 명칭과 문서 수정
- "Grid-point DA" 또는 "1D-Var with ensemble perturbation"
- 전지구 규모에서는 합리적 선택

**옵션 3: 하이브리드**
- 지역 모델: LEnKF
- 전지구 모델: 1D-Var
- 설정 파일로 선택 가능

---

## 7. 참고 문헌

### 코드 파일

1. `lis/core/LIS_lsmMod.F90` (라인 1754-1809)
   - `LIS_lsm_DAmapTileSpaceToObsSpace`

2. `lis/dataassim/algorithm/enkf/enkf_Mod.F90` (라인 365-446)
   - EnKF 메인 루프

3. `lis/dataassim/algorithm/enkf/enkf_general.F90` (라인 217-249)
   - Localization 코드

### 이론

1. Evensen, G. (2003). The Ensemble Kalman Filter: theoretical formulation and practical implementation.
2. Houtekamer & Mitchell (2001). A Sequential Ensemble Kalman Filter for Atmospheric Data Assimilation.
3. Hamill et al. (2001). Distance-Dependent Filtering of Background Error Covariance Estimates.

---

**작성일:** 2025-12-30
**작성자:** Claude Code Analysis
**버전:** LISF 7.5
