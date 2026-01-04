get_gdas.o get_gdas.d : get_gdas.F90
get_gdas.o : LDT_metforcingMod.o
get_gdas.o : gdas_forcingMod.o
get_gdas.o : LDT_constantsMod.o
get_gdas.o : LDT_logMod.o
get_gdas.o : LDT_timeMgrMod.o
get_gdas.o : LDT_coreMod.o
