get_cmap.o get_cmap.d : get_cmap.F90
get_cmap.o : LIS_coreMod.o
get_cmap.o : LIS_constantsMod.o
get_cmap.o : LIS_logMod.o
get_cmap.o : LIS_timeMgrMod.o
get_cmap.o : cmap_forcingMod.o
