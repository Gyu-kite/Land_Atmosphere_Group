read_ecmwf.o read_ecmwf.d : read_ecmwf.F90
read_ecmwf.o : ecmwf_forcingMod.o
read_ecmwf.o : LIS_misc.h
read_ecmwf.o : LIS_metforcingMod.o
read_ecmwf.o : LIS_coreMod.o
read_ecmwf.o : LIS_logMod.o
