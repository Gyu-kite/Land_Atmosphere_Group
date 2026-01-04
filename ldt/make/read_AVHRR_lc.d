read_AVHRR_lc.o read_AVHRR_lc.d : read_AVHRR_lc.F90
read_AVHRR_lc.o : LDT_paramTileInputMod.o
read_AVHRR_lc.o : LDT_logMod.o
read_AVHRR_lc.o : LDT_fileIOMod.o
read_AVHRR_lc.o : LDT_coreMod.o
read_AVHRR_lc.o : LDT_misc.h
read_AVHRR_lc.o : LDT_gridmappingMod.o
