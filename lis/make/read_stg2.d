read_stg2.o read_stg2.d : read_stg2.F90
read_stg2.o : LIS_coreMod.o
read_stg2.o : LIS_logMod.o
read_stg2.o : stg2_forcingMod.o
read_stg2.o : LIS_misc.h
