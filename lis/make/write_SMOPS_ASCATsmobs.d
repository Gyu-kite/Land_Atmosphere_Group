write_SMOPS_ASCATsmobs.o write_SMOPS_ASCATsmobs.d : write_SMOPS_ASCATsmobs.F90
write_SMOPS_ASCATsmobs.o : LIS_coreMod.o
write_SMOPS_ASCATsmobs.o : LIS_historyMod.o
write_SMOPS_ASCATsmobs.o : LIS_fileIOMod.o
write_SMOPS_ASCATsmobs.o : LIS_DAobservationsMod.o
write_SMOPS_ASCATsmobs.o : LIS_logMod.o
write_SMOPS_ASCATsmobs.o : LIS_constantsMod.o
