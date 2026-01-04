readcrd_gdas.o readcrd_gdas.d : readcrd_gdas.F90
readcrd_gdas.o : LDT_logMod.o
readcrd_gdas.o : gdas_forcingMod.o
readcrd_gdas.o : LDT_coreMod.o
