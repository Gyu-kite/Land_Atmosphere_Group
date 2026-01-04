finalize_merra2.o finalize_merra2.d : finalize_merra2.F90
finalize_merra2.o : LIS_coreMod.o
finalize_merra2.o : merra2_forcingMod.o
