get_WRFout.o get_WRFout.d : get_WRFout.F90
get_WRFout.o : LIS_metforcingMod.o
get_WRFout.o : LIS_timeMgrMod.o
get_WRFout.o : LIS_constantsMod.o
get_WRFout.o : LIS_coreMod.o
get_WRFout.o : WRFout_forcingMod.o
get_WRFout.o : LIS_logMod.o
