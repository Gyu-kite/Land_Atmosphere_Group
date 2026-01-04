readUSAFSI.o readUSAFSI.d : readUSAFSI.F90
readUSAFSI.o : LDT_constantsMod.o
readUSAFSI.o : LDT_NetCDF_inc.h
readUSAFSI.o : LDT_coreMod.o
readUSAFSI.o : LDT_misc.h
readUSAFSI.o : LDT_smap_e_oplMod.o
readUSAFSI.o : LDT_logMod.o
