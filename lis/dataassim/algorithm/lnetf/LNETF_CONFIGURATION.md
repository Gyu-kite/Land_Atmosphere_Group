# LNETF Configuration Guide

## Overview

This document describes the configuration parameters for the Local Nonlinear Ensemble Transform Filter (LNETF) data assimilation algorithm in LIS.

## LIS Configuration File Parameters

### Required Parameters

Add these parameters to your `lis.config` file:

```
Data assimilation algorithm:                         "LNETF"
Number of data assimilation instances:               1
Number of ensembles per tile:                        12
```

### LNETF-Specific Parameters

#### 1. Localization Radius Factor

```
LNETF localization radius factor:                    5.0
```

**Description**: Controls the spatial influence radius for localization.

**Default**: 5.0

**Recommended Range**: 3.0 - 10.0

**Guidelines**:
- Smaller values (3.0-5.0): Stricter localization, better for sparse observation networks
- Medium values (5.0-7.0): Balanced, good starting point for most applications
- Larger values (7.0-10.0): Wider influence, suitable for dense observation coverage

**Physical Meaning**:
The effective localization radius in degrees is computed as:
```
radius = sqrt((dx * factor)^2 + (dy * factor)^2)
```
where `dx` and `dy` are the grid resolution in degrees.

**Example**:
- Grid resolution: 0.25° × 0.25°
- Factor: 5.0
- Effective radius: ~1.77° (~200 km at mid-latitudes)

---

## Internal Algorithm Parameters

The following parameters are defined in `lnetf_Mod.F90` and currently use fixed default values. These do not need to be specified in the configuration file for standard usage.

### 2. Forgetting Factor Type (`type_forget`)

**Current Value**: 0 (hardcoded)

**Options**:
- 0: Inflate forecast ensemble (standard approach)
- 1: Inflate forecast ensemble (observed domains only)
- 2: Inflate analysis ensemble
- 3: Inflate analysis ensemble (observed domains only)

**Status**: Fixed value sufficient for most applications

### 3. Forgetting Factor (`forget`)

**Current Value**: 1.0 (hardcoded, no inflation)

**Range**: 1.0 - 1.1

**Meaning**:
- 1.0 = No covariance inflation
- 1.05 = 5% inflation
- 1.1 = 10% inflation

**Status**: Default of 1.0 recommended. EnKF in LIS also does not expose this parameter.

### 4. Transformation Type (`type_trans`)

**Current Value**: 0 (hardcoded)

**Options**:
- 0: Stochastic transformation using random orthonormal matrix (recommended)
- 1: Deterministic transformation

**Status**: Option 0 is the standard LNETF method from Toedter & Ahrens (2015)

### 5. Weights Inflation Type (`type_winf`)

**Current Value**: 0 (hardcoded, disabled)

**Options**:
- 0: No weights inflation
- 1: Inflate weights when N_eff/N > limit_winf

**Status**: Disabled by default. Can prevent particle degeneracy but rarely needed with proper localization.

### 6. Weights Inflation Limit (`limit_winf`)

**Current Value**: 0.0 (hardcoded, disabled)

**Range**: 0.5 - 0.9

**Meaning**: Threshold for effective sample size ratio (N_eff/N) to trigger inflation

**Status**: Only relevant if `type_winf = 1`

---

## Example Configuration

### Basic LNETF Setup for Soil Moisture Assimilation

```bash
#---------------------DATA ASSIMILATION ----------------------------------
Number of data assimilation instances:               1
Data assimilation algorithm:                         "LNETF"
Data assimilation set:                               "AMSR-E(NASA) soil moisture"
Number of state variables:                           4
Data assimilation output interval for diagnostics:  "1da"
Data assimilation number of observation types:       1
Data assimilation output ensemble spread:            1
Data assimilation output processed observations:     1
Data assimilation output innovations:                1

LNETF localization radius factor:                    5.0

Data assimilation scaling strategy:     "CDF matching"
Data assimilation observation domain file:  ../lis_input.d01.nc

# Ensemble configuration
Number of ensembles per tile:              12

# Perturbation settings (same as EnKF)
Perturbations start mode:                 "coldstart"
Forcing perturbation algorithm:           "GMAO scheme"
Forcing perturbation frequency:           "1hr"
State perturbation algorithm:             "GMAO scheme"
State perturbation frequency:             "3hr"
Observation perturbation algorithm:       "GMAO scheme"
Observation perturbation frequency:       "3hr"
```

---

## When to Tune Parameters

### Localization Factor

**Increase factor (6.0 - 10.0) if:**
- Filter is too conservative (small analysis increments)
- Observation network is dense
- Model grid is coarse (> 0.5°)

**Decrease factor (3.0 - 4.0) if:**
- Filter is too aggressive (unrealistic updates)
- Observation network is sparse
- Model grid is fine (< 0.1°)

### Advanced Tuning (Requires Code Modification)

If you need to enable advanced parameters:

1. **Forgetting Factor**: If ensemble spread collapses over time, try `forget = 1.05`
2. **Weights Inflation**: If observing particle degeneracy, enable `type_winf = 1` with `limit_winf = 0.7`

To enable these, modify `lnetf_setup()` in `lnetf_Mod.F90` to read additional parameters:

```fortran
call ESMF_ConfigGetAttribute(LIS_config, &
     lnetf_struc(n,k)%forget, &
     label="LNETF forgetting factor:", rc=status)
```

---

## Comparison with EnKF

| Parameter | EnKF | LNETF | Notes |
|-----------|------|-------|-------|
| Localization | EnKF localization radius factor | LNETF localization radius factor | Same mechanism |
| Inflation | Not configurable in LIS | `forget` (hardcoded 1.0) | Neither exposed to config |
| Transformation | Perturbed observations | Ensemble transform matrix | Different algorithm |
| Non-Gaussianity | Assumes Gaussian | Handles non-Gaussian | LNETF advantage |

---

## References

- Toedter, J., and B. Ahrens, 2015: A Second-Order Exact Ensemble Square Root Filter for Nonlinear Data Assimilation. *Mon. Wea. Rev.*, **143**, 1347-1367.
- PDAF Documentation: http://pdaf.awi.de/
- LIS User Guide: https://lis.gsfc.nasa.gov/

---

## Version History

- 2024-12-28: Initial version for LNETF implementation in LISF 7.5
