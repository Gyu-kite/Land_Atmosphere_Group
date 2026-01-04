readARMObs.o readARMObs.d : readARMObs.F90
readARMObs.o : LVT_timeMgrMod.o
readARMObs.o : map_utils.o
readARMObs.o : ARM_obsMod.o
readARMObs.o : LVT_constantsMod.o
readARMObs.o : LVT_histDataMod.o
readARMObs.o : LVT_misc.h
readARMObs.o : LVT_logMod.o
readARMObs.o : LVT_coreMod.o
