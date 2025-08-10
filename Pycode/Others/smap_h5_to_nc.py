#!/usr/bin/env python3
"""
SMAP L2 HDF5 to NetCDF Converter with Spatial Gridding (Fixed Version)
Convert SMAP L2 soil moisture data from HDF5 swath format to NetCDF regular grid
Fixed to handle integer variables properly
"""

import h5py
import numpy as np
import netCDF4 as nc
import sys
import os
from datetime import datetime
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
import warnings
warnings.filterwarnings('ignore')

def is_valid_data(data, fill_value=-9999.0):
    """
    Check for valid data points, handling both float and integer data types
    
    Parameters:
    -----------
    data : array
        Data array to check
    fill_value : float or int
        Fill value to check against
        
    Returns:
    --------
    valid_mask : array
        Boolean mask of valid data points
    """
    if np.issubdtype(data.dtype, np.integer):
        # For integer data, only check for fill values
        valid_mask = (data != fill_value) & (data != -9999)
    else:
        # For float data, check for NaN and fill values
        valid_mask = ~np.isnan(data) & (data != fill_value) & (data != -9999.0)
    
    return valid_mask

def create_target_grid(lats, lons, target_resolution=0.25):
    """
    Create target grid based on overall lat/lon extent
    
    Parameters:
    -----------
    lats : array
        All latitude values from SMAP L2 file
    lons : array
        All longitude values from SMAP L2 file
    target_resolution : float
        Target grid resolution in degrees
        
    Returns:
    --------
    grid_lats : array
        Gridded latitude array
    grid_lons : array
        Gridded longitude array
    grid_lat_2d : array
        2D latitude grid
    grid_lon_2d : array
        2D longitude grid
    """
    # Use all coordinate data to determine grid extent
    valid_lats = is_valid_data(lats.astype(np.float64))
    valid_lons = is_valid_data(lons.astype(np.float64))
    coord_valid_mask = valid_lats & valid_lons
    
    if np.sum(coord_valid_mask) == 0:
        raise ValueError("No valid coordinate points found")
    
    valid_lat_coords = lats[coord_valid_mask].astype(np.float64)
    valid_lon_coords = lons[coord_valid_mask].astype(np.float64)
    
    # Create target grid based on coordinate extent
    lat_min, lat_max = np.floor(valid_lat_coords.min()), np.ceil(valid_lat_coords.max())
    lon_min, lon_max = np.floor(valid_lon_coords.min()), np.ceil(valid_lon_coords.max())
    
    # Ensure reasonable bounds
    lat_min = max(lat_min, -90)
    lat_max = min(lat_max, 90)
    lon_min = max(lon_min, -180)
    lon_max = min(lon_max, 180)
    
    grid_lats = np.arange(lat_min, lat_max + target_resolution, target_resolution)
    grid_lons = np.arange(lon_min, lon_max + target_resolution, target_resolution)
    
    grid_lon_2d, grid_lat_2d = np.meshgrid(grid_lons, grid_lats)
    
    print(f"Target grid: {len(grid_lats)} x {len(grid_lons)}")
    print(f"Lat range: {lat_min:.2f} to {lat_max:.2f}")
    print(f"Lon range: {lon_min:.2f} to {lon_max:.2f}")
    
    return grid_lats, grid_lons, grid_lat_2d, grid_lon_2d

def grid_smap_data_to_target(lats, lons, data, grid_lat_2d, grid_lon_2d, 
                           target_resolution=0.25, method='nearest'):
    """
    Grid SMAP L2 swath data to predefined target grid
    
    Parameters:
    -----------
    lats : array
        Latitude values from SMAP L2 file
    lons : array
        Longitude values from SMAP L2 file
    data : array
        Data values to be gridded
    grid_lat_2d : array
        Target 2D latitude grid
    grid_lon_2d : array
        Target 2D longitude grid
    target_resolution : float
        Target grid resolution in degrees
    method : str
        Interpolation method ('nearest', 'linear', 'cubic')
    
    Returns:
    --------
    gridded_data : array
        Gridded data array
    """
    
    # Remove invalid data points using improved function
    valid_lats = is_valid_data(lats.astype(np.float64))
    valid_lons = is_valid_data(lons.astype(np.float64))
    valid_data = is_valid_data(data)
    
    # Combined valid mask
    valid_mask = valid_lats & valid_lons & valid_data
    
    if np.sum(valid_mask) == 0:
        print("Warning: No valid data points found")
        # Return grid filled with NaN
        return np.full(grid_lat_2d.shape, np.nan)
    
    valid_lats = lats[valid_mask].astype(np.float64)
    valid_lons = lons[valid_mask].astype(np.float64)
    valid_data_vals = data[valid_mask]
    
    print(f"Valid data points: {len(valid_data_vals)}")
    
    # Grid the data
    try:
        if method == 'nearest':
            # Use KDTree for nearest neighbor (faster for large datasets)
            tree = cKDTree(np.column_stack([valid_lats, valid_lons]))
            distances, indices = tree.query(
                np.column_stack([grid_lat_2d.ravel(), grid_lon_2d.ravel()]),
                k=1
            )
            
            # Set maximum distance threshold (in degrees)
            max_distance = target_resolution * 2
            gridded_data = np.full(grid_lat_2d.shape, np.nan)
            
            valid_interp = distances < max_distance
            gridded_data.ravel()[valid_interp] = valid_data_vals[indices[valid_interp]]
            
        else:
            # Use scipy.interpolate for linear/cubic interpolation
            # Convert integer data to float for interpolation
            if np.issubdtype(valid_data_vals.dtype, np.integer):
                valid_data_vals = valid_data_vals.astype(np.float64)
                
            gridded_data = griddata(
                (valid_lats, valid_lons), valid_data_vals,
                (grid_lat_2d, grid_lon_2d),
                method=method, fill_value=np.nan
            )
            
    except Exception as e:
        print(f"Error during gridding: {e}")
        return np.full(grid_lat_2d.shape, np.nan)
    
    return gridded_data

def convert_smap_l2_to_netcdf(input_h5_file, output_nc_file, 
                              target_resolution=0.25, interp_method='nearest'):
    """
    Convert SMAP L2 HDF5 file to gridded NetCDF format
    
    Parameters:
    -----------
    input_h5_file : str
        Path to input SMAP L2 HDF5 file
    output_nc_file : str
        Path to output NetCDF file
    target_resolution : float
        Target grid resolution in degrees
    interp_method : str
        Interpolation method ('nearest', 'linear', 'cubic')
    """
    
    try:
        # Open HDF5 file
        print(f"Opening HDF5 file: {input_h5_file}")
        with h5py.File(input_h5_file, 'r') as h5_file:
            
            # Print available groups for debugging
            print("Available groups in HDF5 file:")
            for group_name in h5_file.keys():
                print(f"  - {group_name}")
            
            # Access soil moisture retrieval data group
            sm_group_name = "Soil_Moisture_Retrieval_Data"
            if sm_group_name not in h5_file:
                raise ValueError(f"Group '{sm_group_name}' not found in HDF5 file")
            
            sm_group = h5_file[sm_group_name]
            
            # Print available variables for debugging
            print(f"Available variables in {sm_group_name}:")
            for var_name in sm_group.keys():
                var_shape = sm_group[var_name].shape
                var_dtype = sm_group[var_name].dtype
                print(f"  - {var_name}: {var_shape}, {var_dtype}")
            
            # Read coordinate data
            print("Reading coordinate data...")
            if 'latitude' not in sm_group or 'longitude' not in sm_group:
                raise ValueError("Latitude/longitude data not found in HDF5 file")
            
            lats = sm_group['latitude'][:]
            lons = sm_group['longitude'][:]
            
            print(f"Coordinate data shape: {lats.shape}")
            print(f"Data extent: Lat [{np.nanmin(lats):.2f}, {np.nanmax(lats):.2f}], "
                  f"Lon [{np.nanmin(lons):.2f}, {np.nanmax(lons):.2f}]")
            
            # Define variables to process
            variables_to_process = {
                'soil_moisture': {
                    'units': 'cm^3/cm^3',
                    'long_name': 'Soil Moisture',
                    'valid_range': [0.0, 1.0],
                    'dtype': 'f4'
                },
                'soil_moisture_error': {
                    'units': 'cm^3/cm^3',
                    'long_name': 'Soil Moisture Uncertainty',
                    'dtype': 'f4'
                },
                'retrieval_qual_flag': {
                    'units': '1',
                    'long_name': 'Retrieval Quality Flag',
                    'dtype': 'i2'
                },
                'surface_flag': {
                    'units': '1',
                    'long_name': 'Surface Flag',
                    'dtype': 'i2'
                },
                'tb_time_seconds': {
                    'units': 'seconds',
                    'long_name': 'Time in seconds since start of day',
                    'dtype': 'f4'
                }
            }
            
            # Create NetCDF file
            print(f"Creating NetCDF file: {output_nc_file}")
            with nc.Dataset(output_nc_file, 'w', format='NETCDF4') as nc_file:
                
                # Add global attributes
                nc_file.title = "SMAP L2 Soil Moisture Data (Gridded)"
                nc_file.institution = "NASA JPL"
                nc_file.source = "SMAP L2 Radiometer Half-Orbit (Gridded)"
                nc_file.history = f"Converted from swath to grid on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
                nc_file.conventions = "CF-1.6"
                nc_file.grid_resolution = f"{target_resolution} degrees"
                nc_file.interpolation_method = interp_method
                
                # Copy original global attributes
                for attr_name in h5_file.attrs:
                    try:
                        attr_value = h5_file.attrs[attr_name]
                        if isinstance(attr_value, bytes):
                            attr_value = attr_value.decode('utf-8')
                        elif isinstance(attr_value, np.ndarray) and attr_value.dtype.kind == 'S':
                            attr_value = attr_value.astype(str)
                        nc_file.setncattr(f"original_{attr_name}", attr_value)
                    except Exception as e:
                        print(f"Warning: Could not copy attribute {attr_name}: {e}")
                        pass
                
                # Create target grid based on coordinate extent (only once)
                print("Creating target grid based on coordinate extent...")
                grid_lats, grid_lons, grid_lat_2d, grid_lon_2d = create_target_grid(
                    lats, lons, target_resolution
                )
                
                # Create dimensions
                nc_file.createDimension('lat', len(grid_lats))
                nc_file.createDimension('lon', len(grid_lons))
                
                # Create coordinate variables
                lat_var = nc_file.createVariable('lat', 'f4', ('lat',), zlib=True)
                lat_var[:] = grid_lats
                lat_var.units = 'degrees_north'
                lat_var.long_name = 'Latitude'
                lat_var.axis = 'Y'
                
                lon_var = nc_file.createVariable('lon', 'f4', ('lon',), zlib=True)
                lon_var[:] = grid_lons
                lon_var.units = 'degrees_east'
                lon_var.long_name = 'Longitude'
                lon_var.axis = 'X'
                
                # Process each variable
                for var_name, var_attrs in variables_to_process.items():
                    if var_name not in sm_group:
                        print(f"Warning: Variable '{var_name}' not found, skipping...")
                        continue
                    
                    print(f"Processing {var_name}...")
                    
                    # Read data
                    data = sm_group[var_name][:]
                    print(f"  Data shape: {data.shape}, dtype: {data.dtype}")
                    
                    # Grid the data to target grid
                    gridded_data = grid_smap_data_to_target(
                        lats, lons, data, grid_lat_2d, grid_lon_2d, 
                        target_resolution, interp_method
                    )
                    
                    # Determine data type and fill value
                    if var_attrs['dtype'] in ['i2', 'i4']:
                        dtype = var_attrs['dtype']
                        fill_value = -9999
                        # Convert gridded data to integer if needed
                        if not np.issubdtype(gridded_data.dtype, np.integer):
                            # Replace NaN with fill value for integer variables
                            gridded_data = np.where(np.isnan(gridded_data), fill_value, gridded_data)
                            gridded_data = gridded_data.astype(np.int16)
                    else:
                        dtype = 'f4'
                        fill_value = -9999.0
                        # Keep NaN as is for float variables, will be replaced below
                    
                    # Create NetCDF variable
                    nc_var = nc_file.createVariable(
                        var_name, dtype, ('lat', 'lon'), 
                        zlib=True, complevel=6, fill_value=fill_value
                    )
                    
                    # Handle NaN values appropriately
                    if dtype in ['i2', 'i4']:
                        # For integer variables, NaN should already be replaced
                        nc_var[:] = gridded_data
                    else:
                        # For float variables, replace NaN with fill value
                        final_data = np.where(np.isnan(gridded_data), fill_value, gridded_data)
                        nc_var[:] = final_data
                    
                    # Set attributes
                    for attr_name, attr_value in var_attrs.items():
                        if attr_name != 'dtype':  # Skip our internal dtype specification
                            nc_var.setncattr(attr_name, attr_value)
                    
                    # Copy original attributes from HDF5
                    for attr_name in sm_group[var_name].attrs:
                        try:
                            attr_value = sm_group[var_name].attrs[attr_name]
                            if isinstance(attr_value, bytes):
                                attr_value = attr_value.decode('utf-8')
                            elif isinstance(attr_value, np.ndarray) and attr_value.dtype.kind == 'S':
                                attr_value = attr_value.astype(str)
                            nc_var.setncattr(f"original_{attr_name}", attr_value)
                        except Exception as e:
                            print(f"Warning: Could not copy attribute {attr_name} for {var_name}: {e}")
                            pass
                
                # Add data coverage information
                if 'soil_moisture' in [var.name for var in nc_file.variables.values()]:
                    sm_data = nc_file.variables['soil_moisture'][:]
                    if hasattr(sm_data, 'mask'):
                        valid_points = np.sum(~sm_data.mask)
                    else:
                        valid_points = np.sum((sm_data != fill_value) & ~np.isnan(sm_data))
                    total_points = sm_data.size
                    coverage = (valid_points / total_points) * 100
                    
                    nc_file.data_coverage_percent = coverage
                    nc_file.valid_data_points = int(valid_points)
                    nc_file.total_grid_points = int(total_points)
        
        print(f"Successfully converted {input_h5_file} to {output_nc_file}")
        if grid_lats is not None and grid_lons is not None:
            print(f"Grid size: {len(grid_lats)} x {len(grid_lons)}")
            if 'coverage' in locals():
                print(f"Data coverage: {coverage:.1f}%")
        
    except Exception as e:
        print(f"Error during conversion: {str(e)}")
        import traceback
        traceback.print_exc()
        raise

def batch_convert_smap_files(input_dir, output_dir, target_resolution=0.25, 
                           interp_method='nearest', pattern="*.h5"):
    """
    Batch convert multiple SMAP L2 HDF5 files to gridded NetCDF
    """
    import glob
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all HDF5 files
    h5_files = glob.glob(os.path.join(input_dir, pattern))
    
    print(f"Found {len(h5_files)} HDF5 files to convert")
    print(f"Grid resolution: {target_resolution} degrees")
    print(f"Interpolation method: {interp_method}")
    
    successful = 0
    failed = 0
    
    for h5_file in h5_files:
        # Generate output filename
        base_name = os.path.splitext(os.path.basename(h5_file))[0]
        nc_file = os.path.join(output_dir, f"{base_name}_gridded.nc")
        
        try:
            convert_smap_l2_to_netcdf(h5_file, nc_file, target_resolution, interp_method)
            successful += 1
        except Exception as e:
            print(f"Failed to convert {h5_file}: {str(e)}")
            failed += 1
    
    print(f"\nBatch conversion complete:")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")

def main():
    """Main function for command line usage"""
    if len(sys.argv) < 3:
        print("Usage: python smap_h5_to_nc_fixed.py <input_h5_file> <output_nc_file> [resolution] [method]")
        print("   or: python smap_h5_to_nc_fixed.py <input_dir> <output_dir> --batch [resolution] [method]")
        print("")
        print("Options:")
        print("  resolution: Grid resolution in degrees (default: 0.25)")
        print("  method: Interpolation method - 'nearest', 'linear', 'cubic' (default: 'nearest')")
        print("")
        print("Example:")
        print("  python smap_h5_to_nc_fixed.py input.h5 output.nc 0.1 nearest")
        print("  python smap_h5_to_nc_fixed.py /input/dir /output/dir --batch 0.25 linear")
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    
    # Parse optional arguments
    resolution = 0.25
    method = 'nearest'
    
    if '--batch' in sys.argv:
        batch_mode = True
        if len(sys.argv) > 4:
            try:
                resolution = float(sys.argv[4])
            except:
                pass
        if len(sys.argv) > 5:
            method = sys.argv[5]
    else:
        batch_mode = False
        if len(sys.argv) > 3:
            try:
                resolution = float(sys.argv[3])
            except:
                pass
        if len(sys.argv) > 4:
            method = sys.argv[4]
    
    # Validate method
    if method not in ['nearest', 'linear', 'cubic']:
        print(f"Invalid interpolation method: {method}")
        print("Valid methods: 'nearest', 'linear', 'cubic'")
        sys.exit(1)
    
    if batch_mode:
        batch_convert_smap_files(input_path, output_path, resolution, method)
    else:
        convert_smap_l2_to_netcdf(input_path, output_path, resolution, method)

if __name__ == "__main__":
    main()
