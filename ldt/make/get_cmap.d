get_cmap.o get_cmap.d : get_cmap.F90
get_cmap.o : cmap_forcingMod.o
get_cmap.o : LDT_timeMgrMod.o
get_cmap.o : LDT_coreMod.o
get_cmap.o : LDT_constantsMod.o
get_cmap.o : LDT_logMod.o
