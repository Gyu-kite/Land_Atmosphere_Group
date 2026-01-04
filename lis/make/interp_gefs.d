interp_gefs.o interp_gefs.d : interp_gefs.F90
interp_gefs.o : gefs_forcingMod.o
interp_gefs.o : LIS_coreMod.o
interp_gefs.o : LIS_spatialDownscalingMod.o
