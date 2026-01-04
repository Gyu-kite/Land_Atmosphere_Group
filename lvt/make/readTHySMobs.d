readTHySMobs.o readTHySMobs.d : readTHySMobs.F90
readTHySMobs.o : THySM_obsMod.o
readTHySMobs.o : LVT_misc.h
readTHySMobs.o : LVT_logMod.o
readTHySMobs.o : LVT_coreMod.o
readTHySMobs.o : LVT_constantsMod.o
readTHySMobs.o : LVT_timeMgrMod.o
readTHySMobs.o : LVT_histDataMod.o
