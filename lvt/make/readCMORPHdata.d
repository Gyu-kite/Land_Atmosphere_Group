readCMORPHdata.o readCMORPHdata.d : readCMORPHdata.F90
readCMORPHdata.o : LVT_coreMod.o
readCMORPHdata.o : LVT_logMod.o
readCMORPHdata.o : LVT_constantsMod.o
readCMORPHdata.o : LVT_misc.h
readCMORPHdata.o : CMORPH_dataMod.o
readCMORPHdata.o : LVT_timeMgrMod.o
readCMORPHdata.o : LVT_histDataMod.o
