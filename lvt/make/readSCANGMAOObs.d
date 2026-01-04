readSCANGMAOObs.o readSCANGMAOObs.d : readSCANGMAOObs.F90
readSCANGMAOObs.o : SCANGMAO_obsMod.o
readSCANGMAOObs.o : map_utils.o
readSCANGMAOObs.o : LVT_coreMod.o
readSCANGMAOObs.o : LVT_timeMgrMod.o
readSCANGMAOObs.o : LVT_constantsMod.o
readSCANGMAOObs.o : LVT_logMod.o
readSCANGMAOObs.o : LVT_histDataMod.o
