readcrd_ecmwf.o readcrd_ecmwf.d : readcrd_ecmwf.F90
readcrd_ecmwf.o : ecmwf_forcingMod.o
readcrd_ecmwf.o : LIS_logMod.o
readcrd_ecmwf.o : LIS_coreMod.o
