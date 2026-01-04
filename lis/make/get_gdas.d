get_gdas.o get_gdas.d : get_gdas.F90
get_gdas.o : LIS_timeMgrMod.o
get_gdas.o : LIS_constantsMod.o
get_gdas.o : LIS_logMod.o
get_gdas.o : gdas_forcingMod.o
get_gdas.o : LIS_coreMod.o
