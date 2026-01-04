gfs_forcingMod.o gfs_forcingMod.d : gfs_forcingMod.F90
gfs_forcingMod.o : LIS_logMod.o
gfs_forcingMod.o : LIS_timeMgrMod.o
gfs_forcingMod.o : LIS_constantsMod.o
gfs_forcingMod.o : LIS_coreMod.o
