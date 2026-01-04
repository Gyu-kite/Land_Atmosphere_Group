reset_gdas.o reset_gdas.d : reset_gdas.F90
reset_gdas.o : gdas_forcingMod.o
reset_gdas.o : LDT_coreMod.o
reset_gdas.o : LDT_misc.h
reset_gdas.o : LDT_timeMgrMod.o
