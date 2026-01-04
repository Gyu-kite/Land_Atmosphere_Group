read_princeton.o read_princeton.d : read_princeton.F90
read_princeton.o : LIS_forecastMod.o
read_princeton.o : LIS_logMod.o
read_princeton.o : LIS_metforcingMod.o
read_princeton.o : LIS_coreMod.o
read_princeton.o : LIS_misc.h
read_princeton.o : LIS_constantsMod.o
read_princeton.o : LIS_timeMgrMod.o
read_princeton.o : princeton_forcingMod.o
