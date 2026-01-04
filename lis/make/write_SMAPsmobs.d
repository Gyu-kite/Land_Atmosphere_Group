write_SMAPsmobs.o write_SMAPsmobs.d : write_SMAPsmobs.F90
write_SMAPsmobs.o : LIS_fileIOMod.o
write_SMAPsmobs.o : LIS_constantsMod.o
write_SMAPsmobs.o : LIS_historyMod.o
write_SMAPsmobs.o : LIS_coreMod.o
write_SMAPsmobs.o : LIS_logMod.o
