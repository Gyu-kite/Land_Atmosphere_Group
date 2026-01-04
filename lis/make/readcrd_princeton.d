readcrd_princeton.o readcrd_princeton.d : readcrd_princeton.F90
readcrd_princeton.o : LIS_logMod.o
readcrd_princeton.o : LIS_coreMod.o
readcrd_princeton.o : princeton_forcingMod.o
