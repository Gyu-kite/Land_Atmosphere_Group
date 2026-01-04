get_ecmwf.o get_ecmwf.d : get_ecmwf.F90
get_ecmwf.o : LIS_timeMgrMod.o
get_ecmwf.o : LIS_logMod.o
get_ecmwf.o : ecmwf_forcingMod.o
get_ecmwf.o : LIS_coreMod.o
get_ecmwf.o : LIS_constantsMod.o
