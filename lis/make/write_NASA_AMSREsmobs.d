write_NASA_AMSREsmobs.o write_NASA_AMSREsmobs.d : write_NASA_AMSREsmobs.F90
write_NASA_AMSREsmobs.o : LIS_logMod.o
write_NASA_AMSREsmobs.o : LIS_DAobservationsMod.o
write_NASA_AMSREsmobs.o : LIS_constantsMod.o
write_NASA_AMSREsmobs.o : LIS_historyMod.o
write_NASA_AMSREsmobs.o : LIS_fileIOMod.o
write_NASA_AMSREsmobs.o : LIS_coreMod.o
