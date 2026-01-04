readinput_polar.o readinput_polar.d : readinput_polar.F90
readinput_polar.o : map_utils.o
readinput_polar.o : LDT_fileIOMod.o
readinput_polar.o : LDT_logMod.o
readinput_polar.o : LDT_coreMod.o
readinput_polar.o : LDT_domainMod.o
readinput_polar.o : LDT_misc.h
