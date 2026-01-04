read_gldas.o read_gldas.d : read_gldas.F90
read_gldas.o : LIS_coreMod.o
read_gldas.o : LIS_metforcingMod.o
read_gldas.o : LIS_misc.h
read_gldas.o : LIS_logMod.o
read_gldas.o : LIS_timeMgrMod.o
read_gldas.o : gldas_forcingMod.o
