get_princeton.o get_princeton.d : get_princeton.F90
get_princeton.o : princeton_forcingMod.o
get_princeton.o : LIS_coreMod.o
get_princeton.o : LIS_logMod.o
get_princeton.o : LIS_timeMgrMod.o
