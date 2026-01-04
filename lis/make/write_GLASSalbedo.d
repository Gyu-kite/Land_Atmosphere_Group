write_GLASSalbedo.o write_GLASSalbedo.d : write_GLASSalbedo.F90
write_GLASSalbedo.o : LIS_coreMod.o
write_GLASSalbedo.o : LIS_constantsMod.o
write_GLASSalbedo.o : LIS_DAobservationsMod.o
write_GLASSalbedo.o : LIS_logMod.o
write_GLASSalbedo.o : LIS_fileIOMod.o
write_GLASSalbedo.o : LIS_historyMod.o
