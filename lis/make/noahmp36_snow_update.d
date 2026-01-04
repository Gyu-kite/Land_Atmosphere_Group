noahmp36_snow_update.o noahmp36_snow_update.d : noahmp36_snow_update.F90
noahmp36_snow_update.o : module_sf_noahlsm_36.o
noahmp36_snow_update.o : LIS_coreMod.o
noahmp36_snow_update.o : NoahMP36_lsmMod.o
noahmp36_snow_update.o : module_sf_noahmplsm_36.o
