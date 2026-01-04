readGHCNObs.o readGHCNObs.d : readGHCNObs.F90
readGHCNObs.o : LVT_constantsMod.o
readGHCNObs.o : LVT_logMod.o
readGHCNObs.o : GHCN_obsMod.o
readGHCNObs.o : LVT_histDataMod.o
readGHCNObs.o : map_utils.o
readGHCNObs.o : LVT_timeMgrMod.o
readGHCNObs.o : LVT_coreMod.o
