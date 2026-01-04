read_USAFSIobs.o read_USAFSIobs.d : read_USAFSIobs.F90
read_USAFSIobs.o : LIS_timeMgrMod.o
read_USAFSIobs.o : LIS_coreMod.o
read_USAFSIobs.o : LIS_constantsMod.o
read_USAFSIobs.o : LIS_DAobservationsMod.o
read_USAFSIobs.o : LIS_NetCDF_inc.h
read_USAFSIobs.o : LIS_misc.h
read_USAFSIobs.o : LIS_logMod.o
read_USAFSIobs.o : USAFSIobs_Mod.o
