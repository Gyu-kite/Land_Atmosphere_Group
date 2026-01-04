readcrd_nldas20.o readcrd_nldas20.d : readcrd_nldas20.F90
readcrd_nldas20.o : nldas20_forcingMod.o
readcrd_nldas20.o : LIS_logMod.o
readcrd_nldas20.o : LIS_coreMod.o
