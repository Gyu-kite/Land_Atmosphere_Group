read_AVHRR_GFS_lc.o read_AVHRR_GFS_lc.d : read_AVHRR_GFS_lc.F90
read_AVHRR_GFS_lc.o : LDT_logMod.o
read_AVHRR_GFS_lc.o : LDT_gridmappingMod.o
read_AVHRR_GFS_lc.o : LDT_coreMod.o
read_AVHRR_GFS_lc.o : LDT_paramTileInputMod.o
read_AVHRR_GFS_lc.o : LDT_misc.h
read_AVHRR_GFS_lc.o : LDT_fileIOMod.o
