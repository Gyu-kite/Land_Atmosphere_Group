write_SNODEPobs.o write_SNODEPobs.d : write_SNODEPobs.F90
write_SNODEPobs.o : LIS_constantsMod.o
write_SNODEPobs.o : LIS_coreMod.o
write_SNODEPobs.o : LIS_fileIOMod.o
write_SNODEPobs.o : LIS_logMod.o
write_SNODEPobs.o : LIS_DAobservationsMod.o
write_SNODEPobs.o : LIS_historyMod.o
