readWindSatsmObs.o readWindSatsmObs.d : readWindSatsmObs.F90
readWindSatsmObs.o : WindSatsm_obsMod.o
readWindSatsmObs.o : map_utils.o
readWindSatsmObs.o : LDT_constantsMod.o
readWindSatsmObs.o : LDT_coreMod.o
readWindSatsmObs.o : LDT_logMod.o
readWindSatsmObs.o : LDT_timeMgrMod.o
readWindSatsmObs.o : LDT_misc.h
readWindSatsmObs.o : LDT_DAobsDataMod.o
