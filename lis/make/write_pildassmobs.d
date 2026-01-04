write_pildassmobs.o write_pildassmobs.d : write_pildassmobs.F90
write_pildassmobs.o : LIS_coreMod.o
write_pildassmobs.o : LIS_historyMod.o
write_pildassmobs.o : LIS_logMod.o
write_pildassmobs.o : pildassmobs_module.o
write_pildassmobs.o : LIS_DAobservationsMod.o
write_pildassmobs.o : LIS_constantsMod.o
write_pildassmobs.o : LIS_fileIOMod.o
