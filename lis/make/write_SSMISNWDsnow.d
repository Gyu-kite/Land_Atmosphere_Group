write_SSMISNWDsnow.o write_SSMISNWDsnow.d : write_SSMISNWDsnow.F90
write_SSMISNWDsnow.o : LIS_logMod.o
write_SSMISNWDsnow.o : LIS_historyMod.o
write_SSMISNWDsnow.o : LIS_fileIOMod.o
write_SSMISNWDsnow.o : LIS_coreMod.o
write_SSMISNWDsnow.o : LIS_DAobservationsMod.o
write_SSMISNWDsnow.o : LIS_constantsMod.o
