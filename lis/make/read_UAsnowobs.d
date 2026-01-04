read_UAsnowobs.o read_UAsnowobs.d : read_UAsnowobs.F90
read_UAsnowobs.o : LIS_coreMod.o
read_UAsnowobs.o : LIS_misc.h
read_UAsnowobs.o : LIS_logMod.o
read_UAsnowobs.o : UAsnow_obsMod.o
read_UAsnowobs.o : LIS_timeMgrMod.o
read_UAsnowobs.o : LIS_constantsMod.o
read_UAsnowobs.o : map_utils.o
read_UAsnowobs.o : LIS_fileIOMod.o
