readCDFS_GVFObs.o readCDFS_GVFObs.d : readCDFS_GVFObs.F90
readCDFS_GVFObs.o : LDT_misc.h
readCDFS_GVFObs.o : LDT_logMod.o
readCDFS_GVFObs.o : LDT_timeMgrMod.o
readCDFS_GVFObs.o : LDT_coreMod.o
readCDFS_GVFObs.o : CDFS_GVF_obsMod.o
readCDFS_GVFObs.o : LDT_DAobsDataMod.o
readCDFS_GVFObs.o : LDT_constantsMod.o
