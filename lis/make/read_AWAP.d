read_AWAP.o read_AWAP.d : read_AWAP.F90
read_AWAP.o : LIS_logMod.o
read_AWAP.o : LIS_misc.h
read_AWAP.o : AWAP_forcingMod.o
read_AWAP.o : LIS_metforcingMod.o
read_AWAP.o : LIS_coreMod.o
