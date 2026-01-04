read_USGS_lc.o read_USGS_lc.d : read_USGS_lc.F90
read_USGS_lc.o : LDT_logMod.o
read_USGS_lc.o : LDT_paramTileInputMod.o
read_USGS_lc.o : LDT_coreMod.o
read_USGS_lc.o : LDT_fileIOMod.o
read_USGS_lc.o : LDT_misc.h
read_USGS_lc.o : LDT_gridmappingMod.o
