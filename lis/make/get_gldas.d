get_gldas.o get_gldas.d : get_gldas.F90
get_gldas.o : LIS_logMod.o
get_gldas.o : LIS_coreMod.o
get_gldas.o : gldas_forcingMod.o
get_gldas.o : LIS_timeMgrMod.o
