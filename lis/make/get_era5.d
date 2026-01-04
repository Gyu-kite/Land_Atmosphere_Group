get_era5.o get_era5.d : get_era5.F90
get_era5.o : LIS_metforcingMod.o
get_era5.o : LIS_constantsMod.o
get_era5.o : era5_forcingMod.o
get_era5.o : LIS_logMod.o
get_era5.o : LIS_coreMod.o
get_era5.o : LIS_forecastMod.o
get_era5.o : LIS_timeMgrMod.o
