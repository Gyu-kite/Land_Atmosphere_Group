readinput_hrap.o readinput_hrap.d : readinput_hrap.F90
readinput_hrap.o : LDT_domainMod.o
readinput_hrap.o : LDT_coreMod.o
readinput_hrap.o : LDT_xmrg_reader.o
readinput_hrap.o : LDT_fileIOMod.o
readinput_hrap.o : LDT_logMod.o
readinput_hrap.o : LDT_misc.h
readinput_hrap.o : map_utils.o
