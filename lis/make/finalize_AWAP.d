finalize_AWAP.o finalize_AWAP.d : finalize_AWAP.F90
finalize_AWAP.o : LIS_coreMod.o
finalize_AWAP.o : AWAP_forcingMod.o
finalize_AWAP.o : LIS_misc.h
