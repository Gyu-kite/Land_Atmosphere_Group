write_ARSsmobs.o write_ARSsmobs.d : write_ARSsmobs.F90
write_ARSsmobs.o : LIS_coreMod.o
write_ARSsmobs.o : LIS_historyMod.o
write_ARSsmobs.o : LIS_constantsMod.o
write_ARSsmobs.o : LIS_fileIOMod.o
write_ARSsmobs.o : LIS_logMod.o
