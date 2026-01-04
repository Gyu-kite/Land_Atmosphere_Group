read_stg2.o read_stg2.d : read_stg2.F90
read_stg2.o : LDT_logMod.o
read_stg2.o : LDT_misc.h
read_stg2.o : LDT_coreMod.o
read_stg2.o : LDT_metforcingMod.o
read_stg2.o : stg2_forcingMod.o
