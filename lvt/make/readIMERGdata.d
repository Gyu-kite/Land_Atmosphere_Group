readIMERGdata.o readIMERGdata.d : readIMERGdata.F90
readIMERGdata.o : LVT_misc.h
readIMERGdata.o : LVT_constantsMod.o
readIMERGdata.o : IMERG_dataMod.o
readIMERGdata.o : LVT_logMod.o
readIMERGdata.o : LVT_histDataMod.o
readIMERGdata.o : LVT_timeMgrMod.o
readIMERGdata.o : LVT_coreMod.o
