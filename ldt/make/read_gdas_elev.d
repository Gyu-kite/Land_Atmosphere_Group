read_gdas_elev.o read_gdas_elev.d : read_gdas_elev.F90
read_gdas_elev.o : LDT_logMod.o
read_gdas_elev.o : LDT_fileIOMod.o
read_gdas_elev.o : LDT_misc.h
read_gdas_elev.o : gdas_forcingMod.o
read_gdas_elev.o : LDT_coreMod.o
read_gdas_elev.o : LDT_metforcingMod.o
