timeinterp_scan.o timeinterp_scan.d : timeinterp_scan.F90
timeinterp_scan.o : scan_forcingMod.o
timeinterp_scan.o : LIS_logMod.o
timeinterp_scan.o : LIS_coreMod.o
timeinterp_scan.o : LIS_FORC_AttributesMod.o
timeinterp_scan.o : LIS_metforcingMod.o
