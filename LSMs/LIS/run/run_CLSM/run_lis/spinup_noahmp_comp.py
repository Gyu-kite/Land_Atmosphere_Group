#!/usr/bin/env python3
"""
Noah-MP Spin-up Equilibrium Analysis Tool

This script analyzes the last 3 years of a spin-up run to determine if the model
has reached equilibrium. It compares key state variables across years to check
for convergence.

Usage:
    python spinup_equilibrium_check.py <spinup_output_dir>
    
Example:
    python spinup_equilibrium_check.py /path/to/OUTPUT_era5/SURFACEMODEL/spin
"""

import os
import sys
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path

def find_preserved_years(spinup_dir):
    """Find the preserved year directories (last_3, last_2, last_1)"""
    preserved_dirs = {}
    for i in range(1, 4):
        year_dir = os.path.join(spinup_dir, f'last_{i}')
        if os.path.exists(year_dir):
            preserved_dirs[f'Year N-{3-i}' if i > 1 else 'Year N'] = year_dir
    return preserved_dirs

def load_year_data(year_dir, var_name):
    """Load a specific variable from all monthly files in a year"""
    surfacemodel_dir = os.path.join(year_dir, 'SURFACEMODEL')
    if not os.path.exists(surfacemodel_dir):
        return None
    
    # Find all netCDF files
    nc_files = sorted(glob.glob(os.path.join(surfacemodel_dir, '**/*.nc'), recursive=True))
    
    if not nc_files:
        return None
    
    try:
        # Load all files and concatenate
        ds = xr.open_mfdataset(nc_files, combine='by_coords')
        if var_name in ds:
            return ds[var_name]
        else:
            print(f"Warning: Variable {var_name} not found in files")
            return None
    except Exception as e:
        print(f"Error loading data from {year_dir}: {e}")
        return None

def calculate_annual_mean(data):
    """Calculate annual mean, handling potential time dimension issues"""
    if data is None:
        return None
    try:
        if 'time' in data.dims:
            return data.mean(dim='time')
        else:
            return data.mean()
    except Exception as e:
        print(f"Error calculating annual mean: {e}")
        return None

def compare_years(preserved_dirs, var_name):
    """Compare a variable across preserved years"""
    print(f"\nAnalyzing {var_name}...")
    
    annual_means = {}
    for year_label, year_dir in preserved_dirs.items():
        data = load_year_data(year_dir, var_name)
        annual_mean = calculate_annual_mean(data)
        if annual_mean is not None:
            annual_means[year_label] = annual_mean
    
    if len(annual_means) < 2:
        print(f"  Insufficient data to compare {var_name}")
        return None
    
    # Calculate relative differences between consecutive years
    year_labels = sorted(annual_means.keys())
    differences = {}
    
    for i in range(len(year_labels) - 1):
        year1 = year_labels[i]
        year2 = year_labels[i + 1]
        
        mean1 = annual_means[year1]
        mean2 = annual_means[year2]
        
        # Calculate spatial mean of absolute difference
        abs_diff = np.abs(mean2 - mean1)
        spatial_mean_diff = float(abs_diff.mean().values)
        
        # Calculate relative difference (percentage)
        mean_value = float(mean1.mean().values)
        if mean_value != 0:
            rel_diff = (spatial_mean_diff / abs(mean_value)) * 100
        else:
            rel_diff = 0
        
        differences[f"{year1} → {year2}"] = {
            'absolute': spatial_mean_diff,
            'relative': rel_diff,
            'mean_value': mean_value
        }
    
    return differences

def main():
    if len(sys.argv) < 2:
        print("Usage: python spinup_equilibrium_check.py <spinup_output_dir>")
        print("Example: python spinup_equilibrium_check.py /path/to/OUTPUT_era5/SURFACEMODEL/spin")
        sys.exit(1)
    
    spinup_dir = sys.argv[1]
    
    if not os.path.exists(spinup_dir):
        print(f"Error: Directory not found: {spinup_dir}")
        sys.exit(1)
    
    print("=" * 80)
    print("Noah-MP Spin-up Equilibrium Analysis")
    print("=" * 80)
    
    # Find preserved year directories
    preserved_dirs = find_preserved_years(spinup_dir)
    
    if not preserved_dirs:
        print("Error: No preserved year directories found (last_1, last_2, last_3)")
        print("Make sure the spin-up script completed successfully.")
        sys.exit(1)
    
    print(f"\nFound {len(preserved_dirs)} preserved years:")
    for year_label, year_dir in preserved_dirs.items():
        print(f"  {year_label}: {year_dir}")
    
    # Key variables to check for equilibrium
    variables_to_check = [
        'SoilMoist_tavg',   # Soil moisture
        'SoilTemp_tavg',    # Soil temperature
        'LEAFC_tavg',       # Leaf carbon
        'ROOTC_tavg',       # Root carbon
        'LAI_tavg',         # Leaf area index
        'Evap_tavg',        # Evapotranspiration
        'GPP_tavg',         # Gross primary production
    ]
    
    print("\n" + "=" * 80)
    print("Equilibrium Analysis Results")
    print("=" * 80)
    
    all_results = {}
    for var in variables_to_check:
        result = compare_years(preserved_dirs, var)
        if result:
            all_results[var] = result
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    equilibrium_threshold = 1.0  # 1% relative change threshold
    
    for var, differences in all_results.items():
        print(f"\n{var}:")
        is_equilibrium = True
        for transition, values in differences.items():
            status = "✓ EQUILIBRIUM" if values['relative'] < equilibrium_threshold else "✗ NOT EQUILIBRIUM"
            print(f"  {transition}:")
            print(f"    Relative change: {values['relative']:.2f}% {status}")
            print(f"    Absolute change: {values['absolute']:.6f}")
            print(f"    Mean value: {values['mean_value']:.6f}")
            if values['relative'] >= equilibrium_threshold:
                is_equilibrium = False
    
    print("\n" + "=" * 80)
    print("RECOMMENDATION")
    print("=" * 80)
    
    # Check if any variable shows non-equilibrium
    any_non_equilibrium = False
    for var, differences in all_results.items():
        for transition, values in differences.items():
            if values['relative'] >= equilibrium_threshold:
                any_non_equilibrium = True
                break
    
    if any_non_equilibrium:
        print("\n⚠️  MODEL HAS NOT REACHED EQUILIBRIUM")
        print("\nRecommendations:")
        print("1. Run additional spin-up cycles")
        print("2. Check if forcing data is appropriate")
        print("3. Review initial conditions")
        print("\nTo run more spin-up cycles:")
        print("  ./spinup_noahmp.sh [higher_number]")
    else:
        print("\n✓ MODEL HAS REACHED EQUILIBRIUM")
        print("\nYou can proceed with your production run using the final restart file.")
    
    print("\n" + "=" * 80)

if __name__ == '__main__':
    main()
