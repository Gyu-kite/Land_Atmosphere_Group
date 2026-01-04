readcrd_gfs.o readcrd_gfs.d : readcrd_gfs.F90
readcrd_gfs.o : gfs_forcingMod.o
readcrd_gfs.o : LIS_coreMod.o
readcrd_gfs.o : LIS_logMod.o
