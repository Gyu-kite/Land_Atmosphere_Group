finalize_gdas.o finalize_gdas.d : finalize_gdas.F90
finalize_gdas.o : gdas_forcingMod.o
finalize_gdas.o : LIS_coreMod.o
