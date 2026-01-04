readGLDAS1obs.o readGLDAS1obs.d : readGLDAS1obs.F90
readGLDAS1obs.o : GLDAS1obsMod.o
readGLDAS1obs.o : LVT_misc.h
readGLDAS1obs.o : LVT_constantsMod.o
readGLDAS1obs.o : LVT_logMod.o
readGLDAS1obs.o : LVT_coreMod.o
readGLDAS1obs.o : LVT_histDataMod.o
