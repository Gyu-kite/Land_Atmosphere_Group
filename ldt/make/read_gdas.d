read_gdas.o read_gdas.d : read_gdas.F90
read_gdas.o : gdas_forcingMod.o
read_gdas.o : LDT_coreMod.o
read_gdas.o : LDT_constantsMod.o
read_gdas.o : LDT_logMod.o
read_gdas.o : LDT_timeMgrMod.o
read_gdas.o : LDT_misc.h
read_gdas.o : LDT_metforcingMod.o
