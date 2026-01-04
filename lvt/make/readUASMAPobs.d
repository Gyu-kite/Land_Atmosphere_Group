readUASMAPobs.o readUASMAPobs.d : readUASMAPobs.F90
readUASMAPobs.o : LVT_misc.h
readUASMAPobs.o : LVT_coreMod.o
readUASMAPobs.o : LVT_constantsMod.o
readUASMAPobs.o : LVT_timeMgrMod.o
readUASMAPobs.o : LVT_histDataMod.o
readUASMAPobs.o : LVT_logMod.o
readUASMAPobs.o : UASMAP_obsMod.o
