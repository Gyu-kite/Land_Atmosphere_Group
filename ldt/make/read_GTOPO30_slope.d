read_GTOPO30_slope.o read_GTOPO30_slope.d : read_GTOPO30_slope.F90
read_GTOPO30_slope.o : LDT_gridmappingMod.o
read_GTOPO30_slope.o : LDT_logMod.o
read_GTOPO30_slope.o : LDT_fileIOMod.o
read_GTOPO30_slope.o : LDT_paramTileInputMod.o
read_GTOPO30_slope.o : LDT_coreMod.o
read_GTOPO30_slope.o : LDT_misc.h
