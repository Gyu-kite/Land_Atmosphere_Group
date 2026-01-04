get_gfs.o get_gfs.d : get_gfs.F90
get_gfs.o : LDT_metforcingMod.o
get_gfs.o : LDT_logMod.o
get_gfs.o : LDT_coreMod.o
get_gfs.o : LDT_constantsMod.o
get_gfs.o : LDT_timeMgrMod.o
get_gfs.o : gfs_forcingMod.o
