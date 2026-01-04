read_GTOPO30_curv.o read_GTOPO30_curv.d : read_GTOPO30_curv.F90
read_GTOPO30_curv.o : LDT_coreMod.o
read_GTOPO30_curv.o : LDT_misc.h
read_GTOPO30_curv.o : LDT_fileIOMod.o
read_GTOPO30_curv.o : LDT_logMod.o
