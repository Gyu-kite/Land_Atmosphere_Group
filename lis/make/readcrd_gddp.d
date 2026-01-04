readcrd_gddp.o readcrd_gddp.d : readcrd_gddp.F90
readcrd_gddp.o : gddp_forcingMod.o
readcrd_gddp.o : LIS_coreMod.o
readcrd_gddp.o : LIS_logMod.o
