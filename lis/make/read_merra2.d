read_merra2.o read_merra2.d : read_merra2.F90
read_merra2.o : LIS_FORC_AttributesMod.o
read_merra2.o : LIS_coreMod.o
read_merra2.o : LIS_spatialDownscalingMod.o
read_merra2.o : merra2_forcingMod.o
read_merra2.o : LIS_metforcingMod.o
read_merra2.o : LIS_logMod.o
read_merra2.o : LIS_misc.h
