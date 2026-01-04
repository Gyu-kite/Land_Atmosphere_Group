readAmerifluxObs.o readAmerifluxObs.d : readAmerifluxObs.F90
readAmerifluxObs.o : map_utils.o
readAmerifluxObs.o : LVT_logMod.o
readAmerifluxObs.o : LVT_coreMod.o
readAmerifluxObs.o : LVT_histDataMod.o
readAmerifluxObs.o : LVT_constantsMod.o
readAmerifluxObs.o : LVT_timeMgrMod.o
readAmerifluxObs.o : Ameriflux_obsMod.o
