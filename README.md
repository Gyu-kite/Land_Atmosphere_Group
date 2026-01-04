# Multi-Source Data Assimilation in LISF

## Overview

This branch (`multi_source_da`) implements advanced multi-source data assimilation capabilities for the Land Information System Framework (LISF), specifically focusing on **km-based localization** for the Local Nonlinear Ensemble Transform Filter (LNETF).

## Key Features

### 1. LNETF with Km-Based Localization

**Reference Papers:**
- Seo et al. (2021) - JULES + LETKF with σ=30km, patch=150km
- Tak et al. (2025) - Multi-sensor SMAP+ASCAT LETKF
- Toedter and Ahrens (2015) - Second-Order Exact Ensemble Square Root Filter

**Major Improvements:**

#### Km-Based Localization
- **Grid-resolution independent**: Uses actual geographic distances (Haversine)
- **Configurable scale**: Default σ=30km, compact support at 3σ (90km cutoff)
- **Gaspari-Cohn function**: 5th-order polynomial for smooth localization

#### Dual Localization Methods
```fortran
! Method 1: Factor-based (EnKF-compatible)
xcompact = dx * localization_factor
ycompact = dy * localization_factor

! Method 2: Km-based (LETKF paper style) ⭐
dist_km = haversine_km(lon1, lat1, lon2, lat2)
compact_km = 3.0 * localization_scale_km  ! 90km cutoff
weight = gaspari_cohn(dist_km / compact_km)
```

### 2. Multi-Source Simultaneous Assimilation

**Species-Based Framework:**
- Independent observation sources with block-diagonal error covariance
- Zero cross-correlation between different species
- Supports SMAP, ASCAT, and other satellite soil moisture retrievals

**Implementation:**
```fortran
! Observation type with species identifier
type :: obs_type
  integer :: species       ! Unique ID (1=SMAP, 2=ASCAT, etc.)
  real    :: lon, lat
  real    :: value
  real    :: std
  logical :: assim

! Block-diagonal R matrix
if (obs_param(i)%species == obs_param(j)%species) then
  ! Within species: spatial correlation
  cov_factor = obs_param(i)%cross_corr * exp(-distance / corr_len)
else
  ! Different species: independent (cov_factor = 0.0)
endif
```

### 3. Comprehensive Localization Diagnostics

**New Diagnostic Variables:**
```fortran
! Per-grid localization statistics
n_local_obs(:)        ! Count of obs within radius
mean_loc_weight(:)     ! Average weight (0=edge, 1=center)
max_loc_weight(:)      ! Maximum weight (closest obs)
min_loc_weight(:)      ! Minimum weight (farthest obs)
```

### 4. Numerical Stability Enhancements

- **Log-sum-exp trick**: Prevents underflow in likelihood computation
- **NaN checks**: Eigenvalue decomposition safety
- **Effective sample size**: N_eff = 1 / Σ(w²) computation and output

## Algorithm Comparison

### EnKF vs LNETF

| Feature | EnKF | LNETF (This Branch) |
|---------|-------|---------------------|
| Localization | Factor-based only | **Factor + Km-based** |
| Grid-resolution dependent | ❌ Yes | ✅ **No** |
| Academic paper alignment | Standard | ✅ **Seo/Tak papers** |
| Localization diagnostics | ❌ None | ✅ **Full diagnostics** |
| Multi-source support | ✅ Species-based | ✅ **Species-based + simultaneous** |

## Modified Files

### LNETF Implementation
1. `lnetf_Mod.F90` (3,806 lines)
   - Main module with multi-source support
   - Localization diagnostics
   - Km-based and factor-based methods

2. `lnetf_general.F90` (1,280 lines)
   - Core analysis with km-based localization
   - Log-sum-exp likelihood
   - Haversine distance calculation

3. `lnetf_types.F90` (104 lines)
   - Extended type definitions
   - Localization parameters

### EnKF Reference Implementation
4. `enkf_Mod.F90` (1,869 lines)
   - Localization factor support
   - Species-based multi-source

5. `enkf_general.F90` (773 lines)
   - Gaspari-Cohn localization function

6. `enkf_types.F90` (104 lines)
   - Species type definitions

7. `my_lu_decomp.F90` (155 lines)
   - LU decomposition utility

8. `my_matrix_functions.F90` (373 lines)
   - Matrix operations for variance

9. `nr_sort.F90` (102 lines)
   - Numerical Recipes sorting

### Core/Plugins Modifications
10. `LIS_dataAssimMod.F90` (2,125 lines)
    - Data assimilation main logic
    - Multi-source infrastructure

11. `LIS_DAobservationsMod.F90` (2,210 lines)
    - Observation management with species support

12. `LIS_lsmMod.F90` (1,885 lines)
    - Land surface model interface

13. `LIS_dataassim_pluginMod.F90` (238 lines)
    - DA plugin interface

## Configuration

### Km-Based Localization (New)
```fortran
LNETF localization scale (km): 30.0  ! σ = 30km
Compact support: 3*σ = 90km
Weight function: Gaspari-Cohn (5th-order)
```

### Factor-Based Localization (EnKF-compatible)
```fortran
EnKF localization radius factor: 5.0  ! Multiplies dx, dy
```

### Multi-Source Setup
```fortran
! Each observation source has unique species
Species 1: SMAP (σ=30km, R=0.04 m³/m³)
Species 2: ASCAT (σ=45km, R=0.05 m³/m³)

! Block-diagonal error covariance assumed
Cross-correlation (different species): 0.0
```

## Usage

### Compilation
```bash
# Configure LIS
cd lis/
./configure

# Compile with parallel build
./compile -j 8

# Or generate dependencies first (for Cray)
./compile -d
```

### Configuration File
```fortran
! Enable LNETF algorithm
DA algorithm: LNETF

! Multi-source setup
DA observation source: SMAP_L2 ASCAT
DA observation error: 0.04 0.05

! Localization (choose one)
LNETF localization scale (km): 30.0  ! ← Km-based (preferred)
! OR
LNETF localization radius factor: 5.0  ! ← Factor-based (EnKF-style)
```

### Running
```bash
# Navigate to testcase
cd lis/testcases/dataassim/enkf_sm_noah/

# Edit lis.config to point to your data
# Configure multi-source observation paths

# Run LIS
../../LIS lis.config
```

## Diagnostics Output

### Localization Diagnostics (New)
- `n_local_obs`: Number of observations within localization radius
  - Higher values = more observation influence
  - Zero = no observations in radius

- `mean_loc_weight`: Average localization weight
  - Range: 0.0 (edge of radius) to 1.0 (center)
  - Higher = stronger observation influence

- `max_loc_weight`: Maximum weight (closest observation)
  - Identifies most influential observation

- `min_loc_weight`: Minimum weight (farthest observation)
  - Identifies least influential observation

### Effective Sample Size
- `N_eff`: Effective number of ensemble members
  - Calculation: N_eff = 1 / Σ(w_i²)
  - Higher N_eff = better particle diversity
  - N_eff/N_ens ratio indicates filter degeneracy

## References

### Academic Papers
1. **Toedter, J. and Ahrens, B. (2015)**: "A Second-Order Exact Ensemble Square Root Filter for Nonlinear Data Assimilation", *Mon. Wea. Rev.* 143, 1347-1367.
2. **Seo, E., et al. (2021)**: Assimilation of SMAP and ASCAT soil moisture retrievals into JULES land surface model using Local Ensemble Transform Kalman Filter
3. **Tak, Y.-J., et al. (2025)**: Multi-sensor SMAP+ASCAT LETKF implementation
4. **Gaspari, G. and Cohn, S.E. (1999)**: Construction of correlation functions in two and three dimensions, *Q.J.R. Meteorol. Soc.* 125, 523-540.

### LISF Documentation
- Main LISF README: `README.adoc` (top-level)
- LISF User Guide: `docs/LIS_users_guide/`
- Algorithm Documentation: `lis/dataassim/algorithm/CLAUDE.md`
- AGENTS.md: AI coding agent guidelines

## Technical Details

### Localization Weights
```fortran
! Km-based weight calculation
function gaspari_cohn(d)
  ! d = normalized distance (0 to 1)
  ! Compact support: weight = 0.0 for d >= 2.0
  
  if (d <= 1e-3) then
    y = 1.0
  else if (d <= 1.0) then
    ! y = -0.25*d^5 + 0.5*d^4 + 0.625*d^3 - 1.666...*d^2 + 1.0
  else
    ! y = d^5/12.0 - 0.5*d^4 + 0.625*d^3 + 5.0/3.0*d^2 - 5.0*d + 4.0 - 2.0/3.0/d
  endif
  
  gaspari_cohn = y
end function

! Apply to distance
d_norm = dist_km / compact_km  ! compact_km = 3.0 * sigma_km
weight = gaspari_cohn(d_norm)
```

### Effective Sample Size
```fortran
! N_eff = 1 / Σ(w_i²)
! Indicates particle diversity and filter health
! High N_eff/N_ens ratio (> 0.7) is good
! Low ratio (< 0.3) indicates degeneracy
```

## Testing

### Test Cases
- Single source (SMAP only)
- Single source (ASCAT only)
- Multi-source (SMAP + ASCAT simultaneous)
- Varying localization parameters (σ=15km, 30km, 60km)

### Validation Metrics
- Analysis innovation statistics
- Ensemble spread
- RMSE against independent observations
- Localization weight statistics
- Effective sample size evolution

## License

Apache License, Version 2.0

See `LICENSES/` subdirectories for third-party component licenses.

## Authors

- LISF Development Team
- NASA Goddard Space Flight Center

## Changelog

### Version 1.0 (2026-01-04)

**Added:**
- Km-based localization for LNETF (grid-resolution independent)
- Multi-source simultaneous assimilation support
- Comprehensive localization diagnostics (n_local_obs, mean/max/min weights)
- Log-sum-exp trick for numerical stability
- Effective sample size (N_eff) calculation

**Modified:**
- LNETF core algorithms to support both localization methods
- Core/Plugins for multi-source infrastructure

**Acknowledgments:**
- PDAF (Parallel Data Assimilation Framework) - LNETF algorithm design inspiration
- Toedter & Ahrens (2015) - Second-order exact square root filter theory
- Seo et al. (2021) & Tak et al. (2025) - LETKF implementation with km-based localization

## Support

For questions or issues:
- GitHub Discussions: https://github.com/NASA-LIS/LISF/discussions
- LISF Website: https://lis.gsfc.nasa.gov

## Future Work

- [ ] Add multi-source to EnSRF
- [ ] Add km-based localization to EnKF
- [ ] Implement adaptive localization (dynamic σ)
- [ ] Add cross-species correlation estimation
- [ ] Compare simultaneous vs sequential multi-source assimilation
