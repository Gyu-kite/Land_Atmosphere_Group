write_NASASMAPvodobs.o write_NASASMAPvodobs.d : write_NASASMAPvodobs.F90
write_NASASMAPvodobs.o : LIS_historyMod.o
write_NASASMAPvodobs.o : LIS_constantsMod.o
write_NASASMAPvodobs.o : LIS_logMod.o
write_NASASMAPvodobs.o : LIS_fileIOMod.o
write_NASASMAPvodobs.o : LIS_DAobservationsMod.o
write_NASASMAPvodobs.o : LIS_coreMod.o
