get_era5.o get_era5.d : get_era5.F90
get_era5.o : LDT_timeMgrMod.o
get_era5.o : LDT_metforcingMod.o
get_era5.o : era5_forcingMod.o
get_era5.o : LDT_coreMod.o
get_era5.o : LDT_constantsMod.o
get_era5.o : LDT_logMod.o
