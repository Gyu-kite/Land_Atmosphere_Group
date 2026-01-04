readinput_lambert.o readinput_lambert.d : readinput_lambert.F90
readinput_lambert.o : LDT_logMod.o
readinput_lambert.o : LDT_domainMod.o
readinput_lambert.o : LDT_fileIOMod.o
readinput_lambert.o : map_utils.o
readinput_lambert.o : LDT_coreMod.o
readinput_lambert.o : LDT_misc.h
