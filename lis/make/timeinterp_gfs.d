timeinterp_gfs.o timeinterp_gfs.d : timeinterp_gfs.F90
timeinterp_gfs.o : LIS_coreMod.o
timeinterp_gfs.o : LIS_logMod.o
timeinterp_gfs.o : LIS_constantsMod.o
timeinterp_gfs.o : LIS_metforcingMod.o
timeinterp_gfs.o : LIS_FORC_AttributesMod.o
timeinterp_gfs.o : gfs_forcingMod.o
timeinterp_gfs.o : LIS_timeMgrMod.o
