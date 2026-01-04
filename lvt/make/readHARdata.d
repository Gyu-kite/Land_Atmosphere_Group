readHARdata.o readHARdata.d : readHARdata.F90
readHARdata.o : LVT_logMod.o
readHARdata.o : LVT_timeMgrMod.o
readHARdata.o : HAR_dataMod.o
readHARdata.o : LVT_misc.h
readHARdata.o : LVT_histDataMod.o
readHARdata.o : LVT_constantsMod.o
readHARdata.o : LVT_coreMod.o
