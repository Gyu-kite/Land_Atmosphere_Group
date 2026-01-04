write_USAFSIobs.o write_USAFSIobs.d : write_USAFSIobs.F90
write_USAFSIobs.o : LIS_coreMod.o
write_USAFSIobs.o : LIS_constantsMod.o
write_USAFSIobs.o : LIS_fileIOMod.o
write_USAFSIobs.o : LIS_DAobservationsMod.o
write_USAFSIobs.o : LIS_misc.h
write_USAFSIobs.o : LIS_logMod.o
