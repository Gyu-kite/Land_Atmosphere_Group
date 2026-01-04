get_AWRAL.o get_AWRAL.d : get_AWRAL.F90
get_AWRAL.o : AWRAL_forcingMod.o
get_AWRAL.o : LIS_coreMod.o
get_AWRAL.o : LIS_timeMgrMod.o
get_AWRAL.o : LIS_logMod.o
get_AWRAL.o : LIS_constantsMod.o
