read_chirps2.o read_chirps2.d : read_chirps2.F90
read_chirps2.o : LIS_misc.h
read_chirps2.o : LIS_metforcingMod.o
read_chirps2.o : chirps2_forcingMod.o
read_chirps2.o : LIS_logMod.o
read_chirps2.o : LIS_timeMgrMod.o
read_chirps2.o : LIS_coreMod.o
