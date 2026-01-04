get_AWAP.o get_AWAP.d : get_AWAP.F90
get_AWAP.o : LIS_coreMod.o
get_AWAP.o : LIS_constantsMod.o
get_AWAP.o : AWAP_forcingMod.o
get_AWAP.o : LIS_timeMgrMod.o
get_AWAP.o : LIS_logMod.o
