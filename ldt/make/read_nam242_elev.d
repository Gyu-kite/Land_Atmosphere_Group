read_nam242_elev.o read_nam242_elev.d : read_nam242_elev.F90
read_nam242_elev.o : LDT_logMod.o
read_nam242_elev.o : LDT_coreMod.o
read_nam242_elev.o : LDT_metforcingMod.o
read_nam242_elev.o : LDT_fileIOMod.o
read_nam242_elev.o : LDT_misc.h
read_nam242_elev.o : nam242_forcingMod.o
