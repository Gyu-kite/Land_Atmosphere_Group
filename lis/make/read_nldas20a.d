read_nldas20a.o read_nldas20a.d : read_nldas20a.F90
read_nldas20a.o : LIS_metforcingMod.o
read_nldas20a.o : nldas20_forcingMod.o
read_nldas20a.o : LIS_misc.h
read_nldas20a.o : LIS_coreMod.o
read_nldas20a.o : LIS_logMod.o
read_nldas20a.o : LIS_spatialDownscalingMod.o
