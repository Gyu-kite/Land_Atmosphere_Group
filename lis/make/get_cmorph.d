get_cmorph.o get_cmorph.d : get_cmorph.F90
get_cmorph.o : LIS_logMod.o
get_cmorph.o : LIS_coreMod.o
get_cmorph.o : LIS_timeMgrMod.o
get_cmorph.o : LIS_forecastMod.o
get_cmorph.o : cmorph_forcingMod.o
get_cmorph.o : LIS_constantsMod.o
