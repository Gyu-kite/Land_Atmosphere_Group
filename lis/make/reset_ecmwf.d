reset_ecmwf.o reset_ecmwf.d : reset_ecmwf.F90
reset_ecmwf.o : LIS_misc.h
reset_ecmwf.o : ecmwf_forcingMod.o
reset_ecmwf.o : LIS_timeMgrMod.o
reset_ecmwf.o : LIS_coreMod.o
