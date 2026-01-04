read_gdas.o read_gdas.d : read_gdas.F90
read_gdas.o : gdas_forcingMod.o
read_gdas.o : LIS_constantsMod.o
read_gdas.o : LIS_coreMod.o
read_gdas.o : LIS_misc.h
read_gdas.o : LIS_metforcingMod.o
read_gdas.o : LIS_timeMgrMod.o
read_gdas.o : LIS_logMod.o
read_gdas.o : LIS_surfaceModelDataMod.o
