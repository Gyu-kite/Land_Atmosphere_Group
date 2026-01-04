reset_AWAP.o reset_AWAP.d : reset_AWAP.F90
reset_AWAP.o : AWAP_forcingMod.o
reset_AWAP.o : LIS_coreMod.o
reset_AWAP.o : LIS_misc.h
