read_ISMNsmobs.o read_ISMNsmobs.d : read_ISMNsmobs.F90
read_ISMNsmobs.o : LIS_timeMgrMod.o
read_ISMNsmobs.o : LIS_coreMod.o
read_ISMNsmobs.o : LIS_constantsMod.o
read_ISMNsmobs.o : LIS_logMod.o
read_ISMNsmobs.o : map_utils.o
read_ISMNsmobs.o : LIS_misc.h
read_ISMNsmobs.o : ISMNsm_obsMod.o
read_ISMNsmobs.o : LIS_fileIOMod.o
