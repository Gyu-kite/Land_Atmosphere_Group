readcrd_galwem.o readcrd_galwem.d : readcrd_galwem.F90
readcrd_galwem.o : LIS_logMod.o
readcrd_galwem.o : LIS_coreMod.o
readcrd_galwem.o : galwem_forcingMod.o
