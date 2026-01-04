get_WRFoutv2.o get_WRFoutv2.d : get_WRFoutv2.F90
get_WRFoutv2.o : WRFoutv2_forcingMod.o
get_WRFoutv2.o : LIS_timeMgrMod.o
get_WRFoutv2.o : LIS_logMod.o
get_WRFoutv2.o : LIS_coreMod.o
