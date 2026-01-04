get_gswp2.o get_gswp2.d : get_gswp2.F90
get_gswp2.o : LIS_timeMgrMod.o
get_gswp2.o : LIS_coreMod.o
get_gswp2.o : LIS_misc.h
get_gswp2.o : gswp2_forcingMod.o
get_gswp2.o : LIS_logMod.o
get_gswp2.o : LIS_metforcingMod.o
