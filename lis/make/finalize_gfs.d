finalize_gfs.o finalize_gfs.d : finalize_gfs.F90
finalize_gfs.o : LIS_misc.h
finalize_gfs.o : LIS_coreMod.o
finalize_gfs.o : gfs_forcingMod.o
