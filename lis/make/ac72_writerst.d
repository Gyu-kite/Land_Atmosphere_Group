ac72_writerst.o ac72_writerst.d : ac72_writerst.F90
ac72_writerst.o : LIS_constantsMod.o
ac72_writerst.o : LIS_logMod.o
ac72_writerst.o : ac72_lsmMod.o
ac72_writerst.o : LIS_historyMod.o
ac72_writerst.o : LIS_timeMgrMod.o
ac72_writerst.o : LIS_coreMod.o
ac72_writerst.o : LIS_misc.h
ac72_writerst.o : LIS_fileIOMod.o
