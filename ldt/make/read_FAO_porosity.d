read_FAO_porosity.o read_FAO_porosity.d : read_FAO_porosity.F90
read_FAO_porosity.o : LDT_fileIOMod.o
read_FAO_porosity.o : LDT_coreMod.o
read_FAO_porosity.o : LDT_logMod.o
read_FAO_porosity.o : LDT_misc.h
read_FAO_porosity.o : LDT_constantsMod.o
