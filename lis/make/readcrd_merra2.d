readcrd_merra2.o readcrd_merra2.d : readcrd_merra2.F90
readcrd_merra2.o : LIS_logMod.o
readcrd_merra2.o : LIS_coreMod.o
readcrd_merra2.o : merra2_forcingMod.o
