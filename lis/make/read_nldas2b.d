read_nldas2b.o read_nldas2b.d : read_nldas2b.F90
read_nldas2b.o : LIS_coreMod.o
read_nldas2b.o : LIS_metforcingMod.o
read_nldas2b.o : LIS_logMod.o
read_nldas2b.o : LIS_misc.h
read_nldas2b.o : nldas2_forcingMod.o
