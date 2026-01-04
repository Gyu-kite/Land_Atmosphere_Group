read_SRTM_Native_slope.o read_SRTM_Native_slope.d : read_SRTM_Native_slope.F90
read_SRTM_Native_slope.o : LDT_fileIOMod.o
read_SRTM_Native_slope.o : LDT_gridmappingMod.o
read_SRTM_Native_slope.o : LDT_paramTileInputMod.o
read_SRTM_Native_slope.o : LDT_logMod.o
read_SRTM_Native_slope.o : calc_SlopeAspect_module.o
read_SRTM_Native_slope.o : LDT_constantsMod.o
read_SRTM_Native_slope.o : LDT_coreMod.o
read_SRTM_Native_slope.o : LDT_misc.h
