readcrd_snotel.o readcrd_snotel.d : readcrd_snotel.F90
readcrd_snotel.o : LIS_coreMod.o
readcrd_snotel.o : LIS_logMod.o
readcrd_snotel.o : snotel_forcingMod.o
