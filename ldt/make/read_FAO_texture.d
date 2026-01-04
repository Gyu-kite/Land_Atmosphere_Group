read_FAO_texture.o read_FAO_texture.d : read_FAO_texture.F90
read_FAO_texture.o : LDT_misc.h
read_FAO_texture.o : LDT_logMod.o
read_FAO_texture.o : LDT_fileIOMod.o
read_FAO_texture.o : LDT_coreMod.o
