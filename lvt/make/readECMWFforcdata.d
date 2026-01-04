readECMWFforcdata.o readECMWFforcdata.d : readECMWFforcdata.F90
readECMWFforcdata.o : LVT_histDataMod.o
readECMWFforcdata.o : LVT_coreMod.o
readECMWFforcdata.o : ECMWFforc_dataMod.o
readECMWFforcdata.o : LVT_misc.h
readECMWFforcdata.o : LVT_timeMgrMod.o
readECMWFforcdata.o : LVT_constantsMod.o
readECMWFforcdata.o : LVT_logMod.o
