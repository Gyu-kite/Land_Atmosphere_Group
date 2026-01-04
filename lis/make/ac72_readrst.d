ac72_readrst.o ac72_readrst.d : ac72_readrst.F90
ac72_readrst.o : ac72_lsmMod.o
ac72_readrst.o : LIS_coreMod.o
ac72_readrst.o : LIS_misc.h
ac72_readrst.o : LIS_fileIOMod.o
ac72_readrst.o : LIS_constantsMod.o
ac72_readrst.o : LIS_logMod.o
ac72_readrst.o : LIS_timeMgrMod.o
ac72_readrst.o : LIS_historyMod.o
