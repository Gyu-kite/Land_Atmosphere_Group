read_PMW_snow.o read_PMW_snow.d : read_PMW_snow.F90
read_PMW_snow.o : LIS_timeMgrMod.o
read_PMW_snow.o : PMW_snow_Mod.o
read_PMW_snow.o : LIS_misc.h
read_PMW_snow.o : LIS_constantsMod.o
read_PMW_snow.o : LIS_logMod.o
read_PMW_snow.o : LIS_mpiMod.o
read_PMW_snow.o : LIS_pluginIndices.o
read_PMW_snow.o : LIS_coreMod.o
