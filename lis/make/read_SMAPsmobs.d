read_SMAPsmobs.o read_SMAPsmobs.d : read_SMAPsmobs.F90
read_SMAPsmobs.o : SMAPsm_obsMod.o
read_SMAPsmobs.o : LIS_constantsMod.o
read_SMAPsmobs.o : LIS_coreMod.o
read_SMAPsmobs.o : LIS_fileIOMod.o
read_SMAPsmobs.o : LIS_timeMgrMod.o
read_SMAPsmobs.o : LIS_misc.h
read_SMAPsmobs.o : map_utils.o
read_SMAPsmobs.o : LIS_logMod.o
