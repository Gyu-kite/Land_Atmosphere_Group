write_NASASMAPsmobs.o write_NASASMAPsmobs.d : write_NASASMAPsmobs.F90
write_NASASMAPsmobs.o : LIS_DAobservationsMod.o
write_NASASMAPsmobs.o : LIS_logMod.o
write_NASASMAPsmobs.o : LIS_constantsMod.o
write_NASASMAPsmobs.o : LIS_coreMod.o
write_NASASMAPsmobs.o : LIS_historyMod.o
write_NASASMAPsmobs.o : LIS_fileIOMod.o
