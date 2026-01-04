finalize_gdas.o finalize_gdas.d : finalize_gdas.F90
finalize_gdas.o : LDT_coreMod.o
finalize_gdas.o : gdas_forcingMod.o
