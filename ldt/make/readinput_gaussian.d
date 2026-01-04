readinput_gaussian.o readinput_gaussian.d : readinput_gaussian.F90
readinput_gaussian.o : LDT_logMod.o
readinput_gaussian.o : LDT_domainMod.o
readinput_gaussian.o : LDT_coreMod.o
readinput_gaussian.o : map_utils.o
