read_nldas20b.o read_nldas20b.d : read_nldas20b.F90
read_nldas20b.o : LIS_metforcingMod.o
read_nldas20b.o : LIS_logMod.o
read_nldas20b.o : LIS_coreMod.o
read_nldas20b.o : nldas20_forcingMod.o
read_nldas20b.o : LIS_misc.h
