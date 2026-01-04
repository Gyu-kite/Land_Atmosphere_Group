read_gfs_elev.o read_gfs_elev.d : read_gfs_elev.F90
read_gfs_elev.o : gfs_forcingMod.o
read_gfs_elev.o : LDT_metforcingMod.o
read_gfs_elev.o : LDT_coreMod.o
read_gfs_elev.o : LDT_logMod.o
