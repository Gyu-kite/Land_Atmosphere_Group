write_syntheticsmobs.o write_syntheticsmobs.d : write_syntheticsmobs.F90
write_syntheticsmobs.o : LIS_logMod.o
write_syntheticsmobs.o : LIS_fileIOMod.o
write_syntheticsmobs.o : LIS_historyMod.o
write_syntheticsmobs.o : LIS_constantsMod.o
write_syntheticsmobs.o : LIS_DAobservationsMod.o
write_syntheticsmobs.o : syntheticsmobs_module.o
write_syntheticsmobs.o : LIS_coreMod.o
