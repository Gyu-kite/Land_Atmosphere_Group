get_gefs.o get_gefs.d : get_gefs.F90
get_gefs.o : LIS_metforcingMod.o
get_gefs.o : LIS_coreMod.o
get_gefs.o : LIS_constantsMod.o
get_gefs.o : LIS_timeMgrMod.o
get_gefs.o : gefs_forcingMod.o
get_gefs.o : LIS_logMod.o
