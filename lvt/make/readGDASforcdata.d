readGDASforcdata.o readGDASforcdata.d : readGDASforcdata.F90
readGDASforcdata.o : LVT_timeMgrMod.o
readGDASforcdata.o : LVT_coreMod.o
readGDASforcdata.o : GDASforc_dataMod.o
readGDASforcdata.o : LVT_histDataMod.o
readGDASforcdata.o : LVT_constantsMod.o
readGDASforcdata.o : LVT_misc.h
readGDASforcdata.o : LVT_logMod.o
