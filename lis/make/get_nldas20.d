get_nldas20.o get_nldas20.d : get_nldas20.F90
get_nldas20.o : LIS_logMod.o
get_nldas20.o : nldas20_forcingMod.o
get_nldas20.o : LIS_coreMod.o
get_nldas20.o : LIS_constantsMod.o
get_nldas20.o : LIS_metforcingMod.o
get_nldas20.o : LIS_timeMgrMod.o
