readMODSCAGobs.o readMODSCAGobs.d : readMODSCAGobs.F90
readMODSCAGobs.o : LVT_constantsMod.o
readMODSCAGobs.o : LVT_histDataMod.o
readMODSCAGobs.o : LVT_logMod.o
readMODSCAGobs.o : LVT_coreMod.o
readMODSCAGobs.o : LVT_misc.h
readMODSCAGobs.o : map_utils.o
readMODSCAGobs.o : MODSCAG_obsMod.o
readMODSCAGobs.o : LVT_timeMgrMod.o
