write_CDFSgvfobs.o write_CDFSgvfobs.d : write_CDFSgvfobs.F90
write_CDFSgvfobs.o : LIS_historyMod.o
write_CDFSgvfobs.o : LIS_constantsMod.o
write_CDFSgvfobs.o : LIS_fileIOMod.o
write_CDFSgvfobs.o : LIS_logMod.o
write_CDFSgvfobs.o : LIS_DAobservationsMod.o
write_CDFSgvfobs.o : LIS_coreMod.o
