readSMAPvwcobs.o readSMAPvwcobs.d : readSMAPvwcobs.F90
readSMAPvwcobs.o : SMAP_vwcobsMod.o
readSMAPvwcobs.o : LVT_constantsMod.o
readSMAPvwcobs.o : LVT_coreMod.o
readSMAPvwcobs.o : LVT_misc.h
readSMAPvwcobs.o : LVT_logMod.o
readSMAPvwcobs.o : LVT_histDataMod.o
