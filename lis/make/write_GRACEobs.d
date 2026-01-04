write_GRACEobs.o write_GRACEobs.d : write_GRACEobs.F90
write_GRACEobs.o : LIS_DAobservationsMod.o
write_GRACEobs.o : LIS_coreMod.o
write_GRACEobs.o : LIS_fileIOMod.o
write_GRACEobs.o : LIS_logMod.o
write_GRACEobs.o : LIS_constantsMod.o
write_GRACEobs.o : LIS_historyMod.o
