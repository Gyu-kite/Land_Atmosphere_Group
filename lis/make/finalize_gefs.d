finalize_gefs.o finalize_gefs.d : finalize_gefs.F90
finalize_gefs.o : LIS_coreMod.o
finalize_gefs.o : LIS_misc.h
finalize_gefs.o : gefs_forcingMod.o
