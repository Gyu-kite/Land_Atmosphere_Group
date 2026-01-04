readSMOSREXObs.o readSMOSREXObs.d : readSMOSREXObs.F90
readSMOSREXObs.o : LVT_logMod.o
readSMOSREXObs.o : LVT_coreMod.o
readSMOSREXObs.o : LVT_timeMgrMod.o
readSMOSREXObs.o : SMOSREX_obsMod.o
readSMOSREXObs.o : map_utils.o
readSMOSREXObs.o : LVT_histDataMod.o
