readinput_easev2.o readinput_easev2.d : readinput_easev2.F90
readinput_easev2.o : LDT_fileIOMod.o
readinput_easev2.o : LDT_domainMod.o
readinput_easev2.o : LDT_coreMod.o
readinput_easev2.o : LDT_misc.h
readinput_easev2.o : easeV2_utils.o
readinput_easev2.o : LDT_logMod.o
readinput_easev2.o : map_utils.o
