readGLDAS2obs.o readGLDAS2obs.d : readGLDAS2obs.F90
readGLDAS2obs.o : GLDAS2obsMod.o
readGLDAS2obs.o : LVT_coreMod.o
readGLDAS2obs.o : LVT_misc.h
readGLDAS2obs.o : LVT_histDataMod.o
readGLDAS2obs.o : LVT_constantsMod.o
readGLDAS2obs.o : LVT_logMod.o
