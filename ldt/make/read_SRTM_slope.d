read_SRTM_slope.o read_SRTM_slope.d : read_SRTM_slope.F90
read_SRTM_slope.o : LDT_fileIOMod.o
read_SRTM_slope.o : LDT_coreMod.o
read_SRTM_slope.o : LDT_paramTileInputMod.o
read_SRTM_slope.o : LDT_gridmappingMod.o
read_SRTM_slope.o : LDT_logMod.o
read_SRTM_slope.o : LDT_misc.h
