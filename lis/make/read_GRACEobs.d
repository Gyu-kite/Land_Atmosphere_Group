read_GRACEobs.o read_GRACEobs.d : read_GRACEobs.F90
read_GRACEobs.o : LIS_timeMgrMod.o
read_GRACEobs.o : LIS_coreMod.o
read_GRACEobs.o : GRACEobs_module.o
read_GRACEobs.o : LIS_DAobservationsMod.o
read_GRACEobs.o : LIS_constantsMod.o
read_GRACEobs.o : LIS_logMod.o
read_GRACEobs.o : LIS_mpiMod.o
read_GRACEobs.o : LIS_fileIOMod.o
