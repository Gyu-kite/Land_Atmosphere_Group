read_CDFSgvf.o read_CDFSgvf.d : read_CDFSgvf.F90
read_CDFSgvf.o : LIS_logMod.o
read_CDFSgvf.o : CDFSgvf_Mod.o
read_CDFSgvf.o : LIS_timeMgrMod.o
read_CDFSgvf.o : LIS_misc.h
read_CDFSgvf.o : LIS_DAobservationsMod.o
read_CDFSgvf.o : LIS_pluginIndices.o
read_CDFSgvf.o : map_utils.o
read_CDFSgvf.o : LIS_constantsMod.o
read_CDFSgvf.o : LIS_dataAssimMod.o
read_CDFSgvf.o : LIS_mpiMod.o
read_CDFSgvf.o : LIS_coreMod.o
