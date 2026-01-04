read_gddp.o read_gddp.d : read_gddp.F90
read_gddp.o : gddp_forcingMod.o
read_gddp.o : LIS_metforcingMod.o
read_gddp.o : LIS_misc.h
read_gddp.o : LIS_spatialDownscalingMod.o
read_gddp.o : LIS_logMod.o
read_gddp.o : LIS_coreMod.o
read_gddp.o : LIS_constantsMod.o
