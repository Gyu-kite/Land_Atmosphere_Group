write_ISMNsmobs.o write_ISMNsmobs.d : write_ISMNsmobs.F90
write_ISMNsmobs.o : LIS_constantsMod.o
write_ISMNsmobs.o : LIS_historyMod.o
write_ISMNsmobs.o : LIS_fileIOMod.o
write_ISMNsmobs.o : LIS_logMod.o
write_ISMNsmobs.o : LIS_coreMod.o
