timeinterp_gfs.o timeinterp_gfs.d : timeinterp_gfs.F90
timeinterp_gfs.o : LDT_logMod.o
timeinterp_gfs.o : LDT_FORC_AttributesMod.o
timeinterp_gfs.o : LDT_coreMod.o
timeinterp_gfs.o : LDT_timeMgrMod.o
timeinterp_gfs.o : gfs_forcingMod.o
timeinterp_gfs.o : LDT_metforcingMod.o
timeinterp_gfs.o : LDT_constantsMod.o
