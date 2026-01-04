read_ecmwf.o read_ecmwf.d : read_ecmwf.F90
read_ecmwf.o : ecmwf_forcingMod.o
read_ecmwf.o : LDT_metforcingMod.o
read_ecmwf.o : LDT_misc.h
read_ecmwf.o : LDT_coreMod.o
read_ecmwf.o : LDT_logMod.o
