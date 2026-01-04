readASOSWEobs.o readASOSWEobs.d : readASOSWEobs.F90
readASOSWEobs.o : LVT_logMod.o
readASOSWEobs.o : LVT_histDataMod.o
readASOSWEobs.o : LVT_coreMod.o
readASOSWEobs.o : LVT_constantsMod.o
readASOSWEobs.o : ASOSWE_obsMod.o
readASOSWEobs.o : map_utils.o
readASOSWEobs.o : LVT_timeMgrMod.o
readASOSWEobs.o : LVT_misc.h
readASOSWEobs.o : UTM_utils.o
