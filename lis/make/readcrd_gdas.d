readcrd_gdas.o readcrd_gdas.d : readcrd_gdas.F90
readcrd_gdas.o : LIS_logMod.o
readcrd_gdas.o : LIS_coreMod.o
readcrd_gdas.o : gdas_forcingMod.o
