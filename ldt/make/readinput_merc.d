readinput_merc.o readinput_merc.d : readinput_merc.F90
readinput_merc.o : map_utils.o
readinput_merc.o : LDT_fileIOMod.o
readinput_merc.o : LDT_misc.h
readinput_merc.o : LDT_coreMod.o
readinput_merc.o : LDT_logMod.o
readinput_merc.o : LDT_domainMod.o
