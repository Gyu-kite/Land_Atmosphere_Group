readUSGSSFObs.o readUSGSSFObs.d : readUSGSSFObs.F90
readUSGSSFObs.o : LVT_histDataMod.o
readUSGSSFObs.o : LVT_coreMod.o
readUSGSSFObs.o : map_utils.o
readUSGSSFObs.o : USGSSF_obsMod.o
readUSGSSFObs.o : LVT_constantsMod.o
readUSGSSFObs.o : LVT_timeMgrMod.o
readUSGSSFObs.o : LVT_logMod.o
