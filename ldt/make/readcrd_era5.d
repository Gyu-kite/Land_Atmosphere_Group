readcrd_era5.o readcrd_era5.d : readcrd_era5.F90
readcrd_era5.o : era5_forcingMod.o
readcrd_era5.o : LDT_coreMod.o
readcrd_era5.o : LDT_logMod.o
