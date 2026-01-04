readGEOSTEFFObs.o readGEOSTEFFObs.d : readGEOSTEFFObs.F90
readGEOSTEFFObs.o : GEOSTEFF_obsMod.o
readGEOSTEFFObs.o : LDT_constantsMod.o
readGEOSTEFFObs.o : LDT_coreMod.o
readGEOSTEFFObs.o : LDT_logMod.o
readGEOSTEFFObs.o : LDT_DAobsDataMod.o
readGEOSTEFFObs.o : LDT_misc.h
readGEOSTEFFObs.o : LDT_timeMgrMod.o
