read_MODIS_lc.o read_MODIS_lc.d : read_MODIS_lc.F90
read_MODIS_lc.o : LDT_misc.h
read_MODIS_lc.o : LDT_gridmappingMod.o
read_MODIS_lc.o : LDT_logMod.o
read_MODIS_lc.o : LDT_fileIOMod.o
read_MODIS_lc.o : LDT_coreMod.o
read_MODIS_lc.o : LDT_paramTileInputMod.o
