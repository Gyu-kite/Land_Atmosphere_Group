get_gfs.o get_gfs.d : get_gfs.F90
get_gfs.o : LIS_constantsMod.o
get_gfs.o : LIS_logMod.o
get_gfs.o : LIS_timeMgrMod.o
get_gfs.o : gfs_forcingMod.o
get_gfs.o : LIS_coreMod.o
