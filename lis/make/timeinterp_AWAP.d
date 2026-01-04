timeinterp_AWAP.o timeinterp_AWAP.d : timeinterp_AWAP.F90
timeinterp_AWAP.o : LIS_coreMod.o
timeinterp_AWAP.o : LIS_FORC_AttributesMod.o
timeinterp_AWAP.o : LIS_metforcingMod.o
timeinterp_AWAP.o : AWAP_forcingMod.o
timeinterp_AWAP.o : LIS_logMod.o
