read_nam242.o read_nam242.d : read_nam242.F90
read_nam242.o : LIS_misc.h
read_nam242.o : LIS_constantsMod.o
read_nam242.o : LIS_coreMod.o
read_nam242.o : LIS_metforcingMod.o
read_nam242.o : nam242_forcingMod.o
read_nam242.o : LIS_logMod.o
