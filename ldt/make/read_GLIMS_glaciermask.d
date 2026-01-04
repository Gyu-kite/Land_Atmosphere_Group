read_GLIMS_glaciermask.o read_GLIMS_glaciermask.d : read_GLIMS_glaciermask.F90
read_GLIMS_glaciermask.o : LDT_fileIOMod.o
read_GLIMS_glaciermask.o : LDT_misc.h
read_GLIMS_glaciermask.o : LDT_glacierMod.o
read_GLIMS_glaciermask.o : LDT_coreMod.o
read_GLIMS_glaciermask.o : LDT_gridmappingMod.o
read_GLIMS_glaciermask.o : LDT_logMod.o
