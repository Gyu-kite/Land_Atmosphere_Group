readcrd_imerg.o readcrd_imerg.d : readcrd_imerg.F90
readcrd_imerg.o : imerg_forcingMod.o
readcrd_imerg.o : LIS_coreMod.o
readcrd_imerg.o : LIS_logMod.o
