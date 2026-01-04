# AGENTS.md

This file provides guidance for AI coding agents working with the NASA GSFC Land Information System Framework (LISF) Version 7.5.

## Project Overview

LISF is a high-performance land surface modeling framework written primarily in **Fortran 90**. This repository focuses on data assimilation algorithms in `lis/dataassim/algorithm/` that integrate observations into land surface models.

## Build Commands

### Prerequisites
Set required environment variables (see `lis/setup_lisf_env_intel21_parallel.sh` for reference):
```bash
export LIS_ARCH=linux_ifc  # or linux_gfortran
export LIS_FC=mpif90       # Fortran compiler
export LIS_CC=mpicc        # C compiler
# Plus library paths: NetCDF, ESMF, GRIB_API, etc.
```

### Build Process (from project root)
```bash
cd lis/
./configure              # Interactive configuration (generates make/configure.lis)
./compile                # Standard build
./compile -j 8           # Parallel build with 8 jobs
./compile -d             # Generate all dependencies first (recommended for Cray)
./compile -h             # Show all options
```

**Output**: Executable `LIS` in the top-level directory.

### Debugging
```bash
cd lis/make/
make debug               # Print build environment (FFLAGS, VPATH, OBJS, etc.)
make clean               # Clean build artifacts
make realclean           # Full clean (includes dependencies)
```

**Debug flags**: Edit `arch/configure.lis.*` files to add `-g`, `-traceback`, `-check all` (Intel) or `-fbacktrace`, `-fbounds-check` (GNU).

## Testing

**No automated unit tests.** LISF uses a manual **testcase** system:

```bash
cd lis/testcases/dataassim/enkf_sm_noah/  # Example testcase
cat README                                 # Read instructions
# Follow README to download data, edit lis.config, run LIS
../../LIS lis.config                       # Run the test
# Verify outputs manually or compare with expected results
```

**Testcase locations**: `lis/testcases/`, `ldt/testcases/`, `lvt/testcases/`

## Code Style Guidelines

### File Structure
All Fortran modules must include:
```fortran
!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module mymodule_Mod
```

### Module Organization
Algorithms follow a **3-file pattern**:
- `{algorithm}_Mod.F90` - Main module with 6-function interface
- `{algorithm}_types.F90` - Type definitions (observations, parameters)
- `{algorithm}_general.F90` - Core numerical analysis routines

### The 6-Function Interface (MANDATORY for all DA algorithms)
Every algorithm in `lis/dataassim/algorithm/` must implement:
```fortran
subroutine {alg}_init()               ! Allocate structures
subroutine {alg}_setup(k)             ! Initialize from config
subroutine {alg}_increments(n,k)      ! Compute analysis increments
subroutine {alg}_update(n,k)          ! Apply increments to state
subroutine {alg}_diagnostics(n,k)     ! Write output diagnostics
subroutine {alg}_final()              ! Cleanup
```
**Arguments**: `n` = nest index, `k` = DA system index

### Naming Conventions
- **Modules**: `{name}_Mod` (capital M, e.g., `enkf_Mod`, `directInsertion_Mod`)
- **Types**: `{algorithm}_dec` for state containers (e.g., `enkf_dec`, `lnetf_dec`)
- **Local variables**: `snake_case` (e.g., `data_status`, `apply_hadamard`, `n_eff`)
- **Arguments**: `PascalCase` or `Title_Case` (e.g., `N_state`, `Observations`, `State_incr`)
- **Type members**: `snake_case` (e.g., `cat_id`, `region_num`, `fileOpen`)

### Indentation and Formatting
- **Indentation**: 3-5 spaces (NO TABS)
- **Line continuation**: Use `&` at end of line
- **Loop indices**: `i`, `j`, `k` for spatial loops; `n_e` for ensemble members
- **Visibility**: Default `PRIVATE`, explicitly list `public ::` members

### Documentation (BOP/EOP System)
**Every module and subroutine must have a prologue:**
```fortran
!BOP
! !MODULE: enkf_Mod
!
! !DESCRIPTION:
!   Detailed description of the module.
!
! !REVISION HISTORY:
!   27 Feb 2005: Sujay Kumar; Initial Specification
!
! !USES:
!   use ESMF
!   use LIS_coreMod
!EOP
```

For subroutines:
```fortran
!BOP
! !ROUTINE: enkf_increments
! \label{enkf_increments}
!
! !INTERFACE:
subroutine enkf_increments(n,k)
! !DESCRIPTION:
!   Computes analysis increments using EnKF.
!EOP
```

### Preprocessor Directives
- **Required**: `#include "LIS_misc.h"` at top of every module
- **NetCDF conditional**: `#if (defined USE_NETCDF3 || defined USE_NETCDF4)`

### Use Statements
- Group at module level under `! !USES:`
- Use `only` clause for clarity: `use LIS_constantsMod, only : LIS_CONST_PATH_LEN`
- Local `use` statements inside subroutines are acceptable

### Data Structures
- Allocate global structures in `_init()` with dimensions `(LIS_rc%nnest, LIS_rc%ndas)`
- Support multiple nests and multiple simultaneous DA systems
- Example:
  ```fortran
  allocate(enkf_struc(LIS_rc%nnest, LIS_rc%ndas))
  ```

### Error Handling
- Log errors: `write(LIS_logunit,*) '[ERR] Error message'`
- Fatal errors: `call LIS_endrun()`
- Status checks: `call LIS_verify(status, 'Context message')`
- **No silent failures**

### Type Safety
- Always use `implicit none`
- Declare `intent(in)`, `intent(out)`, `intent(inout)` for all arguments
- NO implicit typing

## Common Patterns

### Multi-Source Data Assimilation
EnKF/EnSRF/LNETF support multiple observation sources via `species` identifiers:
```fortran
type :: obs_type
  integer :: species    ! Unique identifier for obs type (e.g., 1=SMAP, 2=ASCAT)
  real    :: lon, lat
  real    :: value
  real    :: std
  logical :: assim
end type
```

### Matrix Operations
Use shared utility `my_matrix_functions.F90`:
- `call row_variance(M, N, A, var)` - Ensemble variance
- NOT `matrix_qr` - use `qr_decomp()` instead

## Anti-Patterns (NEVER DO THIS)

❌ **Suppress type errors**: No `as any` equivalents (this is Fortran, not TypeScript)  
❌ **Skip NASA header**: Every file needs the copyright notice  
❌ **Break the 6-function interface**: All algorithms must implement all 6 functions  
❌ **Forget BOP/EOP**: Documentation tags are mandatory  
❌ **Use tabs**: Always use spaces (3-5)  
❌ **Allocate without `(nnest, ndas)`**: Must support multiple nests/DA systems  
❌ **Commit without testing**: Run `./compile` from root to verify builds  

## Development Workflow

1. **Read existing code first** - Understand the 6-function interface
2. **Check type definitions** - Review `*_types.F90` for data structures
3. **Follow existing patterns** - This is a disciplined codebase
4. **Build from root** - Always test with `./compile` from LISF root
5. **Maintain consistency** - Interface compatibility is critical
6. **Document thoroughly** - Use BOP/EOP for every new function

## File Extensions
- `.F90` - Free-form Fortran 90 (with preprocessing)
- `.f90` - Free-form Fortran 90 (no preprocessing)
- `.F` / `.f` - Fixed-form Fortran (legacy)

## Key Dependencies
- **ESMF** - Earth System Modeling Framework (state management)
- **NetCDF** - I/O for diagnostics (NetCDF3 or NetCDF4)
- **MPI** - Parallel execution
- **LIS Core**: `LIS_coreMod`, `LIS_logMod`, `LIS_DAobservationsMod`

## Additional Resources
- `lis/dataassim/algorithm/CLAUDE.md` - Detailed algorithm documentation
- `docs/LIS_users_guide/` - Comprehensive user guide
- `lis/testcases/` - Example configurations and use cases
