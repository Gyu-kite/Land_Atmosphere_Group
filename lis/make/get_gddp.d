get_gddp.o get_gddp.d : get_gddp.F90
get_gddp.o : LIS_misc.h
get_gddp.o : LIS_coreMod.o
get_gddp.o : LIS_timeMgrMod.o
get_gddp.o : LIS_constantsMod.o
get_gddp.o : LIS_metforcingMod.o
get_gddp.o : gddp_forcingMod.o
get_gddp.o : LIS_logMod.o
