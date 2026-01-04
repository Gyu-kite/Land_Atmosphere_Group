readcrd_gefs.o readcrd_gefs.d : readcrd_gefs.F90
readcrd_gefs.o : gefs_forcingMod.o
readcrd_gefs.o : LIS_coreMod.o
readcrd_gefs.o : LIS_logMod.o
