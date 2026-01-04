finalize_gldas.o finalize_gldas.d : finalize_gldas.F90
finalize_gldas.o : gldas_forcingMod.o
finalize_gldas.o : LIS_misc.h
finalize_gldas.o : LIS_coreMod.o
