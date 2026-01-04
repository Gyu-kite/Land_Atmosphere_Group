read_princeton_elev.o read_princeton_elev.d : read_princeton_elev.F90
read_princeton_elev.o : LDT_fileIOMod.o
read_princeton_elev.o : LDT_metforcingMod.o
read_princeton_elev.o : princeton_forcingMod.o
read_princeton_elev.o : LDT_logMod.o
read_princeton_elev.o : LDT_coreMod.o
