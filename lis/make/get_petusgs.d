get_petusgs.o get_petusgs.d : get_petusgs.F90
get_petusgs.o : LIS_forecastMod.o
get_petusgs.o : LIS_constantsMod.o
get_petusgs.o : LIS_logMod.o
get_petusgs.o : LIS_coreMod.o
get_petusgs.o : petusgs_forcingMod.o
get_petusgs.o : LIS_timeMgrMod.o
