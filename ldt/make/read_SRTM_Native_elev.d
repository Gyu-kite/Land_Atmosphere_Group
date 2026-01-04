read_SRTM_Native_elev.o read_SRTM_Native_elev.d : read_SRTM_Native_elev.F90
read_SRTM_Native_elev.o : LDT_logMod.o
read_SRTM_Native_elev.o : LDT_fileIOMod.o
read_SRTM_Native_elev.o : LDT_misc.h
read_SRTM_Native_elev.o : LDT_gridmappingMod.o
read_SRTM_Native_elev.o : LDT_constantsMod.o
read_SRTM_Native_elev.o : LDT_coreMod.o
read_SRTM_Native_elev.o : LDT_paramTileInputMod.o
