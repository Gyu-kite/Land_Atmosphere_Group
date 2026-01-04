get_narr.o get_narr.d : get_narr.F90
get_narr.o : LIS_timeMgrMod.o
get_narr.o : LIS_logMod.o
get_narr.o : narr_forcingMod.o
get_narr.o : LIS_constantsMod.o
get_narr.o : LIS_coreMod.o
get_narr.o : LIS_metforcingMod.o
