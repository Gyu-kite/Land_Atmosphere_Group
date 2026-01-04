# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This directory contains data assimilation algorithms for the NASA GSFC Land Information System Framework (LISF) Version 7.5. The codebase is written in Fortran 90 and implements six different data assimilation methods for land surface modeling.

## Algorithm Structure

Each algorithm follows a consistent **6-function interface pattern**:

```fortran
subroutine algorithm_init()               ! Allocate structures
subroutine algorithm_setup(k)             ! Initialize from config
subroutine algorithm_increments(n,k)      ! Compute analysis increments
subroutine algorithm_update(n,k)          ! Apply increments to state
subroutine algorithm_diagnostics(n,k)     ! Write output diagnostics
subroutine algorithm_final()              ! Cleanup
```

## Algorithm Implementations

### 1. Direct Insertion (DI) - `di/`
- **Main Module**: `directInsertion_Mod.F90` (201 lines)
- **Description**: Simplest DA method that directly replaces model state with observed data
- **Use Case**: Basic observation assimilation without statistical filtering

### 2. Ensemble Kalman Filter (EKF) - `ekf/`
- **Main Module**: `ekf_Mod.F90` (2,133 lines)
- **Description**: Screen-level assimilation using ensemble Kalman filtering (Met Office approach)
- **Dependencies**: ESMF, NetCDF

### 3. GMAO Ensemble Kalman Filter (EnKF) - `enkf/`
- **Main Module**: `enkf_Mod.F90` (1,609 lines)
- **Description**: Full ensemble Kalman filter from NASA GMAO (Rolf Reichle)
- **Supporting Modules**:
  - `enkf_types.F90` - Type definitions for observations
  - `enkf_general.F90` - Core `enkf_analysis()` subroutine
  - `my_matrix_functions.F90` - Matrix operations (variance, QR decomposition)
  - `my_lu_decomp.F90` - LU decomposition
  - `nr_sort.F90` - Numerical recipes sorting
- **Key Features**:
  - Multi-source capability via species-based framework
  - Block-diagonal error covariance for independent observation sources
  - Localization support via `localization_factor` parameter

### 4. Ensemble Square Root Filter (EnSRF) - `ensrf/`
- **Main Module**: `ensrf_Mod.F90` (1,784 lines)
- **Description**: Computationally efficient variant of EnKF avoiding explicit error covariance
- **Supporting Modules**:
  - `ensrf_types.F90` - Type definitions
  - `ensrf_general.F90` - Core `ensrf_analysis()` subroutine
  - Shared: `my_matrix_functions.F90`
- **Advantages**: Better numerical stability via square-root filters

### 5. Ensemble Kalman Smoother GRACE (EnKSGRACE) - `enksgrace/`
- **Main Module**: `enksgrace_Mod.F90` (1,470 lines)
- **Description**: EnKF variant optimized for GRACE satellite data
- **Supporting Modules**:
  - `enks_types.F90` - Type definitions
  - `enks_general.F90` - Core `enks_analysis()` subroutine
  - `enks_matrix_functions.F90` - Matrix operations with GRACE-specific handling
  - `enks_lu_decomp.F90` - LU decomposition
  - `nr_fft.F90` - Fast Fourier Transform utilities
- **Special Requirements**: Requires 30-minute windows as multiples of LIS timestep

### 6. Particle Filter (PF) - `pf/`
- **Main Module**: `pf_Mod.F90` (1,563 lines)
- **Description**: Non-linear DA using sequential importance resampling
- **Supporting Modules**:
  - `pf_types.F90` - Type definitions
  - `pf_general.F90` - Core `pf_analysis()` subroutine
  - Shared: `my_matrix_functions.F90`

### 7. Local Nonlinear Ensemble Transform Filter (LNETF) - `lnetf/` ⭐ NEW
- **Main Module**: `lnetf_Mod.F90` (~1,600 lines)
- **Description**: Domain-localized nonlinear ensemble transform filter for non-Gaussian distributions
- **Algorithm Reference**: J. Toedter and B. Ahrens, "A Second-Order Exact Ensemble Square Root Filter for Nonlinear Data Assimilation", *Mon. Wea. Rev.* 143 (2015) 1347
- **Implementation Source**: Adapted from PDAF (Parallel Data Assimilation Framework) by Lars Nerger and Paul Kirchgessner
- **Supporting Modules**:
  - `lnetf_types.F90` - Type definitions (identical structure to enkf_types)
  - `lnetf_general.F90` - Core `lnetf_analysis()` subroutine
  - Shared: `my_matrix_functions.F90`
- **Key Features**:
  - **Particle weights**: Uses likelihood-based weights instead of Kalman gain
  - **Transform matrix**: Computes A = diag(w) - w*w^T from normalized weights
  - **Square-root filter**: Uses eigenvalue decomposition for numerical stability
  - **Nonlinear capability**: Better handles non-Gaussian error distributions than EnKF
  - **Localization**: Supports spatial localization via `localization_factor` parameter
- **LNETF-specific Parameters** (in `lnetf_dec` type):
  ```fortran
  integer :: type_forget = 0       ! Forgetting factor type
  integer :: type_trans = 0        ! Ensemble transformation type
  integer :: type_winf = 0         ! Weights inflation type
  real :: forget = 1.0             ! Forgetting factor (covariance inflation)
  real :: limit_winf = 0.0         ! Limit for weights inflation
  real :: eff_sample_size = 0.0    ! Effective sample size diagnostic
  ```
- **Algorithm Steps**:
  1. Compute particle weights (likelihood for each ensemble member)
  2. Normalize weights and optionally inflate if N_eff/N > limit_winf
  3. Compute transform matrix A = diag(w) - w*w^T
  4. Calculate square root via eigenvalue decomposition
  5. Apply random rotation (or deterministic transformation)
  6. Transform ensemble: X^a = X^f * W
- **When to Use LNETF vs EnKF**:
  - LNETF: Strongly non-Gaussian errors, nonlinear observation operators
  - EnKF: Approximately Gaussian errors, linear or mildly nonlinear systems
- **Status**: ⚠️ Initial implementation adapted from EnKF template. Core PDAF LNETF analysis routine needs integration.

## Core Data Structures

All ensemble/filter algorithms use **identical observation type definitions**:

```fortran
type obs_type
  integer :: species        ! Observation type identifier
  real    :: lon, lat       ! Location
  real    :: value          ! Observation value
  real    :: std            ! Error standard deviation
  logical :: assim          ! Assimilation flag
end type

type obs_param_type
  integer :: species
  character(LEN=*) :: path
  real    :: std            ! Error std dev
  real    :: cross_corr     ! Cross-correlation parameter
  real    :: corr_len       ! Correlation length
end type

type update_region_type
  real :: lat_min, lat_max
  real :: lon_min, lon_max
end type
```

## Multi-Source Data Assimilation

The **EnKF implementation supports multi-source assimilation** via a species-based framework:

- Each observation type gets a unique `species` identifier
- Error covariances are **block-diagonal** across species (independent sources)
- Within-species spatial correlations are supported via `cross_corr` and `corr_len` parameters
- Example: Simultaneous SMAP + ASCAT soil moisture assimilation is supported

**Key code pattern** (from `enkf_Mod.F90:687-820`):
```fortran
! R matrix construction handles multiple species
if (obs_param(i)%species == obs_param(j)%species) then
  ! Within species: apply spatial correlation
  cov_factor = obs_param(i)%cross_corr * exp(-distance / obs_param(i)%corr_len)
else
  ! Different species: block-diagonal (no correlation)
  cov_factor = 0.0
endif
```

## Build System

**No local Makefiles** - algorithms are compiled as part of the parent LIS framework.

### Compilation Process

From the top-level directory (`/land1/user/gychoi/LIS/test_merge_DA_LISF/`):

```bash
# Configure LIS (creates make/Makefile)
./configure

# Compile LIS (builds all algorithms)
./compile

# Advanced options
./compile -j 8              # Parallel build with 8 jobs
./compile -d                # Generate all dependencies first (recommended for Cray)
./compile -h                # Show all options
```

**Build outputs**:
- Executable: `LIS` (in top-level directory)
- Object files: `make/*.o` (including algorithm modules)

### Compilation Order

Dependencies are automatically resolved by the build system:
```
[Algorithm Main Module] → [Algorithm Types] → [Algorithm General]
  → [Shared Utilities] → [LIS Framework] → [ESMF/NetCDF]
```

## Key Dependencies

- **ESMF** (Earth System Modeling Framework) - State/observation management
- **NetCDF** - I/O for diagnostics (when `USE_NETCDF3` or `USE_NETCDF4` defined)
- **LIS Core Modules**:
  - `LIS_coreMod` - Core runtime control (`LIS_rc`)
  - `LIS_DAobservationsMod` - Observation state handling
  - `LIS_surfaceModelMod` - Land surface model interface
  - `LIS_fileIOMod`, `LIS_historyMod`, `LIS_timeMgrMod` - I/O and time management

## Important Coding Conventions

### Preprocessor Directives
- `#include "LIS_misc.h"` - Required at top of all modules
- `#if (defined USE_NETCDF3 || defined USE_NETCDF4)` - NetCDF conditional compilation

### Data Structures
- Multi-dimensional allocations indexed as `(LIS_rc%nnest, LIS_rc%ndas)` for spatial domains and DA systems
- All algorithms support multiple nests and multiple simultaneous DA systems

### Matrix Operations
- Shared utility `my_matrix_functions.F90` provides:
  - `get_variance()` - Ensemble variance calculation
  - `qr_decomp()` - QR decomposition (not matrix_qr)
  - Hadamard product support for covariance localization

### File Naming Convention
- Main module: `<algorithm>_Mod.F90` (capitalized "M")
- Types: `<algorithm>_types.F90`
- Core analysis: `<algorithm>_general.F90`
- EnKF-specific utilities: `my_*.F90` (GMAO convention)
- EnKSGRACE-specific utilities: `enks_*.F90`

## Configuration Files

Algorithm selection and parameters are specified in LIS configuration files (not in this directory). Key parameters include:

- Data assimilation algorithm (selects which algorithm to use)
- Observation sources and their `species` identifiers
- Error statistics per observation type (`std`, `cross_corr`, `corr_len`)
- Localization parameters (`localization_factor`)
- Update regions (geographic bounding boxes)

## Documentation Notes

- See `enkf/LENKF_다중소스_분석.md` for detailed Korean-language analysis of LENKF multi-source assimilation capabilities
- NASA license header required in all files (see existing files for template)
- Revision history tracked in module header comments

## Development Workflow

When modifying algorithms:

1. **Read the algorithm's main module first** to understand the 6-function interface
2. **Check type definitions** in `*_types.F90` to understand data structures
3. **Review the analysis subroutine** in `*_general.F90` for core algorithm logic
4. **Test from parent directory**: Always build and test using `./compile` from LISF root
5. **Maintain interface consistency**: All algorithms must provide the same 6 functions
6. **Multi-nest/multi-DA support**: Ensure all allocations handle `(nnest, ndas)` indexing

## Current Work

Based on recent files, active development areas include:
- **LNETF Implementation** (Dec 2024) - New algorithm in `lnetf/` directory
  - Initial framework created from EnKF template
  - Core modules: `lnetf_Mod.F90`, `lnetf_types.F90`, `lnetf_general.F90`
  - Based on PDAF implementation at `/home/gychoi/PDAF/src`
  - **Next steps**:
    - Integrate PDAF likelihood-based weight calculation (from `PDAF_lnetf_analysis.F90`)
    - Replace Kalman gain update with transform matrix approach
    - Add eigenvalue decomposition for square-root filter
    - Implement weights inflation and effective sample size diagnostics
    - Test with nonlinear observation operators
- Multi-source data assimilation (SMAP + ASCAT) - see `enkf/LENKF_다중소스_분석.md`
- EnKF algorithm improvements - see backup files `enkf_Mod.F90.bak*` (dated Dec 2024)

## PDAF Integration Reference

LNETF implementation references PDAF source files:
- `/home/gychoi/PDAF/src/PDAF_lnetf.F90` - Module parameters and initialization
- `/home/gychoi/PDAF/src/PDAF_lnetf_analysis.F90` - Core analysis routine with weight calculation
- `/home/gychoi/PDAF/src/PDAF_lnetf_update.F90` - Update control and loop structure

Key PDAF LNETF components to integrate:
1. **Likelihood computation** (lines 183-201 in PDAF_lnetf_analysis.F90): Call user-defined `U_likelihood_l()` for each member
2. **Weight normalization** (lines 218-233): Normalize and handle zero-weight cases
3. **Transform matrix construction** (lines 250-260): A = diag(w) - w*w^T
4. **Eigenvalue decomposition** (lines 290-299): `syev` for computing sqrt(A)
5. **Ensemble transformation** (lines 362-367): T = sqrt(m) * sqrt(A) * rndmat + w
