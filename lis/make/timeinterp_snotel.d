timeinterp_snotel.o timeinterp_snotel.d : timeinterp_snotel.F90
timeinterp_snotel.o : LIS_logMod.o
timeinterp_snotel.o : LIS_FORC_AttributesMod.o
timeinterp_snotel.o : snotel_forcingMod.o
timeinterp_snotel.o : LIS_coreMod.o
timeinterp_snotel.o : LIS_metforcingMod.o
