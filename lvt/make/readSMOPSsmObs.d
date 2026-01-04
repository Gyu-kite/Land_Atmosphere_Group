readSMOPSsmObs.o readSMOPSsmObs.d : readSMOPSsmObs.F90
readSMOPSsmObs.o : SMOPSsm_obsMod.o
readSMOPSsmObs.o : LVT_histDataMod.o
readSMOPSsmObs.o : LVT_constantsMod.o
readSMOPSsmObs.o : LVT_coreMod.o
readSMOPSsmObs.o : LVT_logMod.o
readSMOPSsmObs.o : map_utils.o
readSMOPSsmObs.o : LVT_timeMgrMod.o
readSMOPSsmObs.o : LVT_misc.h
