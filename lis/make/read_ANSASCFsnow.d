read_ANSASCFsnow.o read_ANSASCFsnow.d : read_ANSASCFsnow.F90
read_ANSASCFsnow.o : LIS_constantsMod.o
read_ANSASCFsnow.o : ANSASCFsnow_Mod.o
read_ANSASCFsnow.o : LIS_misc.h
read_ANSASCFsnow.o : LIS_mpiMod.o
read_ANSASCFsnow.o : LIS_timeMgrMod.o
read_ANSASCFsnow.o : LIS_logMod.o
read_ANSASCFsnow.o : LIS_coreMod.o
