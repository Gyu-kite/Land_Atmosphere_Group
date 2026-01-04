get_nldas2.o get_nldas2.d : get_nldas2.F90
get_nldas2.o : nldas2_forcingMod.o
get_nldas2.o : LIS_coreMod.o
get_nldas2.o : LIS_logMod.o
get_nldas2.o : LIS_timeMgrMod.o
get_nldas2.o : LIS_metforcingMod.o
get_nldas2.o : LIS_constantsMod.o
