readcrd_ecmwf.o readcrd_ecmwf.d : readcrd_ecmwf.F90
readcrd_ecmwf.o : LDT_coreMod.o
readcrd_ecmwf.o : ecmwf_forcingMod.o
readcrd_ecmwf.o : LDT_logMod.o
