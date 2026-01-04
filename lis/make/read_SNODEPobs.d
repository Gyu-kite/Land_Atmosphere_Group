read_SNODEPobs.o read_SNODEPobs.d : read_SNODEPobs.F90
read_SNODEPobs.o : LIS_timeMgrMod.o
read_SNODEPobs.o : LIS_DAobservationsMod.o
read_SNODEPobs.o : SNODEPobs_Mod.o
read_SNODEPobs.o : LIS_mpiMod.o
read_SNODEPobs.o : LIS_historyMod.o
read_SNODEPobs.o : LIS_misc.h
read_SNODEPobs.o : LIS_constantsMod.o
read_SNODEPobs.o : LIS_logMod.o
read_SNODEPobs.o : LIS_coreMod.o
