get_snotel.o get_snotel.d : get_snotel.F90
get_snotel.o : snotel_forcingMod.o
get_snotel.o : LIS_metforcingMod.o
get_snotel.o : LIS_coreMod.o
get_snotel.o : LIS_logMod.o
get_snotel.o : LIS_timeMgrMod.o
