readLISrunoffdata.o readLISrunoffdata.d : readLISrunoffdata.F90
readLISrunoffdata.o : LIS_constantsMod.o
readLISrunoffdata.o : LIS_fileIOMod.o
readLISrunoffdata.o : LIS_coreMod.o
readLISrunoffdata.o : LISrunoffdataMod.o
readLISrunoffdata.o : LIS_misc.h
readLISrunoffdata.o : LIS_logMod.o
