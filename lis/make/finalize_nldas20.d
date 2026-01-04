finalize_nldas20.o finalize_nldas20.d : finalize_nldas20.F90
finalize_nldas20.o : LIS_coreMod.o
finalize_nldas20.o : nldas20_forcingMod.o
