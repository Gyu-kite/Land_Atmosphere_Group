readGRACEObs.o readGRACEObs.d : readGRACEObs.F90
readGRACEObs.o : GRACE_obsMod.o
readGRACEObs.o : LVT_timeMgrMod.o
readGRACEObs.o : LVT_logMod.o
readGRACEObs.o : LVT_misc.h
readGRACEObs.o : LVT_histDataMod.o
readGRACEObs.o : LVT_coreMod.o
