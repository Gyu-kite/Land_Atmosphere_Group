write_THySM.o write_THySM.d : write_THySM.F90
write_THySM.o : LIS_coreMod.o
write_THySM.o : LIS_logMod.o
write_THySM.o : LIS_constantsMod.o
write_THySM.o : LIS_fileIOMod.o
write_THySM.o : LIS_historyMod.o
write_THySM.o : LIS_DAobservationsMod.o
