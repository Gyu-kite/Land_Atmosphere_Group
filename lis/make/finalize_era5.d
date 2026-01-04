finalize_era5.o finalize_era5.d : finalize_era5.F90
finalize_era5.o : era5_forcingMod.o
finalize_era5.o : LIS_coreMod.o
