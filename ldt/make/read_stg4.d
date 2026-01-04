read_stg4.o read_stg4.d : read_stg4.F90
read_stg4.o : LDT_coreMod.o
read_stg4.o : LDT_metforcingMod.o
read_stg4.o : LDT_misc.h
read_stg4.o : LDT_logMod.o
read_stg4.o : LDT_constantsMod.o
read_stg4.o : stg4_forcingMod.o
