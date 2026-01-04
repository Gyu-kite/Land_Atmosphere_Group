readUSDMobs.o readUSDMobs.d : readUSDMobs.F90
readUSDMobs.o : LVT_histDataMod.o
readUSDMobs.o : map_utils.o
readUSDMobs.o : LVT_misc.h
readUSDMobs.o : LVT_coreMod.o
readUSDMobs.o : LVT_logMod.o
readUSDMobs.o : USDM_obsMod.o
readUSDMobs.o : LVT_constantsMod.o
