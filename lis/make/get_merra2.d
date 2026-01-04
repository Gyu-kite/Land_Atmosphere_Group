get_merra2.o get_merra2.d : get_merra2.F90
get_merra2.o : merra2_forcingMod.o
get_merra2.o : LIS_logMod.o
get_merra2.o : LIS_timeMgrMod.o
get_merra2.o : LIS_metforcingMod.o
get_merra2.o : LIS_coreMod.o
get_merra2.o : LIS_forecastMod.o
get_merra2.o : LIS_constantsMod.o
