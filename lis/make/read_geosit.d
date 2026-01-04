read_geosit.o read_geosit.d : read_geosit.F90
read_geosit.o : LIS_metforcingMod.o
read_geosit.o : geosit_forcingMod.o
read_geosit.o : LIS_misc.h
read_geosit.o : LIS_spatialDownscalingMod.o
read_geosit.o : LIS_logMod.o
read_geosit.o : LIS_coreMod.o
read_geosit.o : LIS_FORC_AttributesMod.o
