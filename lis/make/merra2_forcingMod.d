merra2_forcingMod.o merra2_forcingMod.d : merra2_forcingMod.F90
merra2_forcingMod.o : LIS_constantsMod.o
merra2_forcingMod.o : LIS_spatialDownscalingMod.o
merra2_forcingMod.o : LIS_coreMod.o
merra2_forcingMod.o : LIS_forecastMod.o
merra2_forcingMod.o : LIS_logMod.o
merra2_forcingMod.o : LIS_timeMgrMod.o
