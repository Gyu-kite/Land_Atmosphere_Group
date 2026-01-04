read_ARSsmobs.o read_ARSsmobs.d : read_ARSsmobs.F90
read_ARSsmobs.o : LIS_logMod.o
read_ARSsmobs.o : LIS_fileIOMod.o
read_ARSsmobs.o : LIS_timeMgrMod.o
read_ARSsmobs.o : ARSsm_obsMod.o
read_ARSsmobs.o : LIS_constantsMod.o
read_ARSsmobs.o : LIS_coreMod.o
