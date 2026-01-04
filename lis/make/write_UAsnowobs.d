write_UAsnowobs.o write_UAsnowobs.d : write_UAsnowobs.F90
write_UAsnowobs.o : LIS_logMod.o
write_UAsnowobs.o : LIS_historyMod.o
write_UAsnowobs.o : LIS_constantsMod.o
write_UAsnowobs.o : LIS_coreMod.o
write_UAsnowobs.o : LIS_fileIOMod.o
