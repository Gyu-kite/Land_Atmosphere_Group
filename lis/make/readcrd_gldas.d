readcrd_gldas.o readcrd_gldas.d : readcrd_gldas.F90
readcrd_gldas.o : LIS_logMod.o
readcrd_gldas.o : gldas_forcingMod.o
readcrd_gldas.o : LIS_coreMod.o
