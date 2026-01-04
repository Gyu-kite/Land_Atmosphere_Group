readSNODASObs.o readSNODASObs.d : readSNODASObs.F90
readSNODASObs.o : LVT_histDataMod.o
readSNODASObs.o : LVT_logMod.o
readSNODASObs.o : LVT_timeMgrMod.o
readSNODASObs.o : LVT_coreMod.o
readSNODASObs.o : SNODAS_obsMod.o
readSNODASObs.o : LVT_constantsMod.o
