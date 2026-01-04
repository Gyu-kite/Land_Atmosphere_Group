get_imerg.o get_imerg.d : get_imerg.F90
get_imerg.o : LIS_forecastMod.o
get_imerg.o : LIS_timeMgrMod.o
get_imerg.o : LIS_constantsMod.o
get_imerg.o : imerg_forcingMod.o
get_imerg.o : LIS_coreMod.o
get_imerg.o : LIS_logMod.o
