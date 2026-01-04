write_syntheticwlobs.o write_syntheticwlobs.d : write_syntheticwlobs.F90
write_syntheticwlobs.o : LIS_historyMod.o
write_syntheticwlobs.o : LIS_DAobservationsMod.o
write_syntheticwlobs.o : LIS_constantsMod.o
write_syntheticwlobs.o : LIS_fileIOMod.o
write_syntheticwlobs.o : LIS_logMod.o
write_syntheticwlobs.o : LIS_coreMod.o
write_syntheticwlobs.o : syntheticwlobs_module.o
