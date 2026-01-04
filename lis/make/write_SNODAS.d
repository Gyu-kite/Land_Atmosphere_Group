write_SNODAS.o write_SNODAS.d : write_SNODAS.F90
write_SNODAS.o : LIS_coreMod.o
write_SNODAS.o : LIS_DAobservationsMod.o
write_SNODAS.o : LIS_constantsMod.o
write_SNODAS.o : LIS_fileIOMod.o
write_SNODAS.o : LIS_historyMod.o
write_SNODAS.o : LIS_logMod.o
