readUWETObs.o readUWETObs.d : readUWETObs.F90
readUWETObs.o : LVT_constantsMod.o
readUWETObs.o : LVT_misc.h
readUWETObs.o : LVT_logMod.o
readUWETObs.o : LVT_coreMod.o
readUWETObs.o : UWET_obsMod.o
readUWETObs.o : LVT_histDataMod.o
