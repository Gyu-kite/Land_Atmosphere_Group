read_princeton.o read_princeton.d : read_princeton.F90
read_princeton.o : LDT_constantsMod.o
read_princeton.o : LDT_logMod.o
read_princeton.o : princeton_forcingMod.o
read_princeton.o : LDT_misc.h
read_princeton.o : LDT_timeMgrMod.o
read_princeton.o : LDT_metforcingMod.o
read_princeton.o : LDT_coreMod.o
