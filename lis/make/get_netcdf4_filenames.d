get_netcdf4_filenames.o get_netcdf4_filenames.d : get_netcdf4_filenames.F90
get_netcdf4_filenames.o : LIS_coreMod.o
get_netcdf4_filenames.o : LIS_forecastMod.o
get_netcdf4_filenames.o : LIS_logMod.o
