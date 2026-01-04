readGRDCobs.o readGRDCobs.d : readGRDCobs.F90
readGRDCobs.o : LVT_constantsMod.o
readGRDCobs.o : GRDC_obsMod.o
readGRDCobs.o : LVT_histDataMod.o
readGRDCobs.o : map_utils.o
readGRDCobs.o : LVT_logMod.o
readGRDCobs.o : LVT_coreMod.o
readGRDCobs.o : LVT_timeMgrMod.o
