read_stg4.o read_stg4.d : read_stg4.F90
read_stg4.o : stg4_forcingMod.o
read_stg4.o : LIS_logMod.o
read_stg4.o : LIS_metforcingMod.o
read_stg4.o : LIS_coreMod.o
read_stg4.o : LIS_misc.h
