get_agrradps.o get_agrradps.d : get_agrradps.F90
get_agrradps.o : LIS_timeMgrMod.o
get_agrradps.o : agrradps_forcingMod.o
get_agrradps.o : LIS_coreMod.o
