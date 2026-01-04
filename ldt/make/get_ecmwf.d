get_ecmwf.o get_ecmwf.d : get_ecmwf.F90
get_ecmwf.o : LDT_constantsMod.o
get_ecmwf.o : LDT_metforcingMod.o
get_ecmwf.o : LDT_timeMgrMod.o
get_ecmwf.o : ecmwf_forcingMod.o
get_ecmwf.o : LDT_logMod.o
get_ecmwf.o : LDT_coreMod.o
