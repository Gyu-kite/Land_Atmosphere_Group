gfs_forcingMod.o gfs_forcingMod.d : gfs_forcingMod.F90
gfs_forcingMod.o : LDT_logMod.o
gfs_forcingMod.o : LDT_timeMgrMod.o
gfs_forcingMod.o : LDT_coreMod.o
gfs_forcingMod.o : LDT_constantsMod.o
