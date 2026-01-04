readGLERLobs.o readGLERLobs.d : readGLERLobs.F90
readGLERLobs.o : LVT_timeMgrMod.o
readGLERLobs.o : LVT_constantsMod.o
readGLERLobs.o : LVT_misc.h
readGLERLobs.o : LVT_histDataMod.o
readGLERLobs.o : GLERL_dataMod.o
readGLERLobs.o : LVT_logMod.o
readGLERLobs.o : map_utils.o
readGLERLobs.o : LVT_coreMod.o
