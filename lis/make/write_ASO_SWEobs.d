write_ASO_SWEobs.o write_ASO_SWEobs.d : write_ASO_SWEobs.F90
write_ASO_SWEobs.o : LIS_fileIOMod.o
write_ASO_SWEobs.o : LIS_DAobservationsMod.o
write_ASO_SWEobs.o : LIS_constantsMod.o
write_ASO_SWEobs.o : LIS_historyMod.o
write_ASO_SWEobs.o : LIS_logMod.o
write_ASO_SWEobs.o : LIS_coreMod.o
