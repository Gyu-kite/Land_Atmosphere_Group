timeinterp_petusgs.o timeinterp_petusgs.d : timeinterp_petusgs.F90
timeinterp_petusgs.o : LIS_metforcingMod.o
timeinterp_petusgs.o : LIS_coreMod.o
timeinterp_petusgs.o : petusgs_forcingMod.o
timeinterp_petusgs.o : LIS_FORC_AttributesMod.o
timeinterp_petusgs.o : LIS_logMod.o
timeinterp_petusgs.o : LIS_forecastMod.o
