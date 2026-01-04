read_nldas2a.o read_nldas2a.d : read_nldas2a.F90
read_nldas2a.o : LIS_spatialDownscalingMod.o
read_nldas2a.o : LIS_misc.h
read_nldas2a.o : LIS_metforcingMod.o
read_nldas2a.o : LIS_coreMod.o
read_nldas2a.o : nldas2_forcingMod.o
read_nldas2a.o : LIS_logMod.o
