readGHCNANNdata.o readGHCNANNdata.d : readGHCNANNdata.F90
readGHCNANNdata.o : map_utils.o
readGHCNANNdata.o : LDT_logMod.o
readGHCNANNdata.o : LDT_misc.h
readGHCNANNdata.o : GHCN_ANNdataMod.o
readGHCNANNdata.o : LDT_ANNMod.o
readGHCNANNdata.o : LDT_timeMgrMod.o
readGHCNANNdata.o : LDT_constantsMod.o
readGHCNANNdata.o : LDT_coreMod.o
