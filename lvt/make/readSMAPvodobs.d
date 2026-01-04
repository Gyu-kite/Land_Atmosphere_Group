readSMAPvodobs.o readSMAPvodobs.d : readSMAPvodobs.F90
readSMAPvodobs.o : LVT_histDataMod.o
readSMAPvodobs.o : LVT_coreMod.o
readSMAPvodobs.o : LVT_misc.h
readSMAPvodobs.o : LVT_timeMgrMod.o
readSMAPvodobs.o : LVT_constantsMod.o
readSMAPvodobs.o : LVT_logMod.o
readSMAPvodobs.o : SMAP_vodobsMod.o
