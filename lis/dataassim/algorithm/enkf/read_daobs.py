#!/usr/bin/env python3
"""
Read LIS DAOBS files (.1gs4r format) and convert to NetCDF

Usage:
    python read_daobs.py input.1gs4r output.nc
"""

import numpy as np
from scipy.io import FortranFile
import xarray as xr
import sys
import os

def read_daobs_file(filename, nc=1440, nr=720):
    """
    Read LIS DAOBS Fortran unformatted binary file

    Parameters:
    -----------
    filename : str
        Path to .1gs4r file
    nc : int
        Number of columns (longitude points)
    nr : int
        Number of rows (latitude points)

    Returns:
    --------
    smobs_unsc : 2D array
        Unscaled soil moisture observations
    smobs : 2D array
        Scaled soil moisture observations
    """

    print(f"Reading {filename}...")
    print(f"Grid size: {nc} × {nr}")

    # Open Fortran unformatted file (big-endian)
    f = FortranFile(filename, 'r', header_dtype=np.dtype('>i4'))

    # Read two 2D arrays
    # Note: Fortran arrays are column-major, so we need to transpose
    # Data is also big-endian
    smobs_unsc = f.read_reals(dtype=np.dtype('>f4')).reshape((nc, nr), order='F').T
    smobs = f.read_reals(dtype=np.dtype('>f4')).reshape((nc, nr), order='F').T

    f.close()

    print(f"  smobs_unsc: min={smobs_unsc.min():.4f}, max={smobs_unsc.max():.4f}")
    print(f"  smobs:      min={smobs.min():.4f}, max={smobs.max():.4f}")

    # Count non-fill values (assuming -9999.0 is the fill value)
    valid_unsc = np.sum(smobs_unsc > -9000)
    valid_sc = np.sum(smobs > -9000)
    print(f"  Valid pixels: {valid_unsc} (unscaled), {valid_sc} (scaled)")

    return smobs_unsc, smobs

def create_netcdf(smobs_unsc, smobs, output_file, nc=1440, nr=720):
    """
    Create NetCDF file from DAOBS data

    Parameters:
    -----------
    smobs_unsc : 2D array
        Unscaled observations
    smobs : 2D array
        Scaled observations
    output_file : str
        Output NetCDF filename
    nc : int
        Number of longitude points
    nr : int
        Number of latitude points
    """

    # Create coordinate arrays (0.25 degree global grid)
    # Assuming grid starts at -180, -90
    lon = np.linspace(-180 + 0.125, 180 - 0.125, nc)
    lat = np.linspace(-90 + 0.125, 90 - 0.125, nr)

    # Create xarray Dataset
    ds = xr.Dataset({
        'smobs_unscaled': (['lat', 'lon'], smobs_unsc),
        'smobs_scaled': (['lat', 'lon'], smobs)
    },
    coords={
        'lon': lon,
        'lat': lat
    })

    # Add attributes
    ds['smobs_unscaled'].attrs = {
        'long_name': 'Unscaled soil moisture observations',
        'units': 'm3/m3',
        '_FillValue': -9999.0
    }

    ds['smobs_scaled'].attrs = {
        'long_name': 'Scaled soil moisture observations',
        'units': 'm3/m3',
        '_FillValue': -9999.0
    }

    ds.attrs = {
        'title': 'LIS Data Assimilation Observations',
        'source': 'LIS DAOBS output',
        'Conventions': 'CF-1.6'
    }

    # Save to NetCDF
    print(f"\nSaving to {output_file}...")
    ds.to_netcdf(output_file)
    print("Done!")

    return ds

def main():
    if len(sys.argv) < 2:
        print("Usage: python read_daobs.py input.1gs4r [output.nc]")
        print("\nExample:")
        print("  python read_daobs.py 202306052330.1gs4r")
        print("  python read_daobs.py 202306052330.1gs4r output.nc")
        sys.exit(1)

    input_file = sys.argv[1]

    if not os.path.exists(input_file):
        print(f"Error: File not found: {input_file}")
        sys.exit(1)

    # Set output filename
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        output_file = os.path.splitext(input_file)[0] + '.nc'

    # Read DAOBS file
    smobs_unsc, smobs = read_daobs_file(input_file)

    # Create NetCDF
    ds = create_netcdf(smobs_unsc, smobs, output_file)

    print(f"\nYou can now view with:")
    print(f"  ncview {output_file}")
    print(f"  ncdump -h {output_file}")

if __name__ == '__main__':
    main()
