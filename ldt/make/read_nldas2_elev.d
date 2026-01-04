read_nldas2_elev.o read_nldas2_elev.d : read_nldas2_elev.F90
read_nldas2_elev.o : LDT_fileIOMod.o
read_nldas2_elev.o : LDT_metforcingMod.o
read_nldas2_elev.o : LDT_coreMod.o
read_nldas2_elev.o : LDT_logMod.o
read_nldas2_elev.o : nldas2_forcingMod.o
