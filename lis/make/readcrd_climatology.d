readcrd_climatology.o readcrd_climatology.d : readcrd_climatology.F90
readcrd_climatology.o : LIS_coreMod.o
readcrd_climatology.o : climatology_forcingMod.o
readcrd_climatology.o : LIS_logMod.o
