readGlobSnowObs.o readGlobSnowObs.d : readGlobSnowObs.F90
readGlobSnowObs.o : LVT_constantsMod.o
readGlobSnowObs.o : LVT_misc.h
readGlobSnowObs.o : LVT_histDataMod.o
readGlobSnowObs.o : GlobSnow_obsMod.o
readGlobSnowObs.o : LVT_logMod.o
readGlobSnowObs.o : LVT_coreMod.o
