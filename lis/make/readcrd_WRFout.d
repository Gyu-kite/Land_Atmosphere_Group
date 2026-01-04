readcrd_WRFout.o readcrd_WRFout.d : readcrd_WRFout.F90
readcrd_WRFout.o : LIS_logMod.o
readcrd_WRFout.o : WRFout_forcingMod.o
readcrd_WRFout.o : LIS_coreMod.o
