read_ANSASNWDsnow.o read_ANSASNWDsnow.d : read_ANSASNWDsnow.F90
read_ANSASNWDsnow.o : LIS_coreMod.o
read_ANSASNWDsnow.o : ANSASNWDsnow_Mod.o
read_ANSASNWDsnow.o : LIS_DAobservationsMod.o
read_ANSASNWDsnow.o : LIS_timeMgrMod.o
read_ANSASNWDsnow.o : LIS_misc.h
read_ANSASNWDsnow.o : LIS_mpiMod.o
read_ANSASNWDsnow.o : LIS_constantsMod.o
read_ANSASNWDsnow.o : LIS_logMod.o
read_ANSASNWDsnow.o : LIS_pluginIndices.o
