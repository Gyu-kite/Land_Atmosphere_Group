readcrd_AWAP.o readcrd_AWAP.d : readcrd_AWAP.F90
readcrd_AWAP.o : AWAP_forcingMod.o
readcrd_AWAP.o : LIS_logMod.o
readcrd_AWAP.o : LIS_coreMod.o
