read_gefs_reforecast.o read_gefs_reforecast.d : read_gefs_reforecast.F90
read_gefs_reforecast.o : LIS_misc.h
read_gefs_reforecast.o : LIS_logMod.o
read_gefs_reforecast.o : LIS_coreMod.o
read_gefs_reforecast.o : LIS_metforcingMod.o
read_gefs_reforecast.o : gefs_forcingMod.o
