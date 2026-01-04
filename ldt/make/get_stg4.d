get_stg4.o get_stg4.d : get_stg4.F90
get_stg4.o : stg4_forcingMod.o
get_stg4.o : LDT_coreMod.o
get_stg4.o : LDT_timeMgrMod.o
get_stg4.o : LDT_logMod.o
get_stg4.o : LDT_constantsMod.o
