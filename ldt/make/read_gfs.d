read_gfs.o read_gfs.d : read_gfs.F90
read_gfs.o : LDT_constantsMod.o
read_gfs.o : LDT_timeMgrMod.o
read_gfs.o : LDT_misc.h
read_gfs.o : gfs_forcingMod.o
read_gfs.o : LDT_coreMod.o
read_gfs.o : LDT_logMod.o
read_gfs.o : LDT_metforcingMod.o
