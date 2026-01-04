readinput_latlon.o readinput_latlon.d : readinput_latlon.F90
readinput_latlon.o : map_utils.o
readinput_latlon.o : LDT_misc.h
readinput_latlon.o : LDT_fileIOMod.o
readinput_latlon.o : LDT_domainMod.o
readinput_latlon.o : LDT_coreMod.o
readinput_latlon.o : LDT_logMod.o
