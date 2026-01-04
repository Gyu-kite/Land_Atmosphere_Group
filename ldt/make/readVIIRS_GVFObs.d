readVIIRS_GVFObs.o readVIIRS_GVFObs.d : readVIIRS_GVFObs.F90
readVIIRS_GVFObs.o : LDT_logMod.o
readVIIRS_GVFObs.o : LDT_misc.h
readVIIRS_GVFObs.o : LDT_coreMod.o
readVIIRS_GVFObs.o : LDT_timeMgrMod.o
readVIIRS_GVFObs.o : LDT_DAobsDataMod.o
readVIIRS_GVFObs.o : VIIRS_GVF_obsMod.o
readVIIRS_GVFObs.o : LDT_constantsMod.o
