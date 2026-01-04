read_SRTM_elev.o read_SRTM_elev.d : read_SRTM_elev.F90
read_SRTM_elev.o : LDT_misc.h
read_SRTM_elev.o : LDT_logMod.o
read_SRTM_elev.o : LDT_coreMod.o
read_SRTM_elev.o : LDT_paramTileInputMod.o
read_SRTM_elev.o : LDT_fileIOMod.o
read_SRTM_elev.o : LDT_gridmappingMod.o
