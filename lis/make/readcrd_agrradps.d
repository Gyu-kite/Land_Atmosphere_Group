readcrd_agrradps.o readcrd_agrradps.d : readcrd_agrradps.F90
readcrd_agrradps.o : agrradps_forcingMod.o
readcrd_agrradps.o : LIS_coreMod.o
readcrd_agrradps.o : LIS_logMod.o
