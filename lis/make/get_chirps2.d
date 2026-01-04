get_chirps2.o get_chirps2.d : get_chirps2.F90
get_chirps2.o : LIS_constantsMod.o
get_chirps2.o : LIS_coreMod.o
get_chirps2.o : LIS_timeMgrMod.o
get_chirps2.o : LIS_logMod.o
get_chirps2.o : LIS_forecastMod.o
get_chirps2.o : chirps2_forcingMod.o
