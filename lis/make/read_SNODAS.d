read_SNODAS.o read_SNODAS.d : read_SNODAS.F90
read_SNODAS.o : LIS_pluginIndices.o
read_SNODAS.o : LIS_timeMgrMod.o
read_SNODAS.o : map_utils.o
read_SNODAS.o : SNODAS_Mod.o
read_SNODAS.o : LIS_DAobservationsMod.o
read_SNODAS.o : LIS_coreMod.o
read_SNODAS.o : LIS_dataAssimMod.o
read_SNODAS.o : LIS_misc.h
read_SNODAS.o : LIS_mpiMod.o
read_SNODAS.o : LIS_logMod.o
read_SNODAS.o : LIS_constantsMod.o
