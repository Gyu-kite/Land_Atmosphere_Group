read_era5.o read_era5.d : read_era5.F90
read_era5.o : LIS_logMod.o
read_era5.o : LIS_FORC_AttributesMod.o
read_era5.o : LIS_coreMod.o
read_era5.o : LIS_metforcingMod.o
read_era5.o : LIS_spatialDownscalingMod.o
read_era5.o : LIS_misc.h
read_era5.o : era5_forcingMod.o
