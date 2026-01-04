readsyntheticsmObs.o readsyntheticsmObs.d : readsyntheticsmObs.F90
readsyntheticsmObs.o : syntheticsm_obsMod.o
readsyntheticsmObs.o : LDT_coreMod.o
readsyntheticsmObs.o : LDT_DAobsDataMod.o
readsyntheticsmObs.o : LDT_misc.h
readsyntheticsmObs.o : LDT_logMod.o
readsyntheticsmObs.o : LDT_constantsMod.o
readsyntheticsmObs.o : map_utils.o
readsyntheticsmObs.o : LDT_timeMgrMod.o
