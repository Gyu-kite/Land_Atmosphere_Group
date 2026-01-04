readcrd_narr.o readcrd_narr.d : readcrd_narr.F90
readcrd_narr.o : LIS_logMod.o
readcrd_narr.o : narr_forcingMod.o
readcrd_narr.o : LIS_coreMod.o
