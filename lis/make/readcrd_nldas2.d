readcrd_nldas2.o readcrd_nldas2.d : readcrd_nldas2.F90
readcrd_nldas2.o : nldas2_forcingMod.o
readcrd_nldas2.o : LIS_logMod.o
readcrd_nldas2.o : LIS_coreMod.o
