write_ANSASCFsnow.o write_ANSASCFsnow.d : write_ANSASCFsnow.F90
write_ANSASCFsnow.o : LIS_coreMod.o
write_ANSASCFsnow.o : LIS_fileIOMod.o
write_ANSASCFsnow.o : LIS_constantsMod.o
write_ANSASCFsnow.o : LIS_logMod.o
write_ANSASCFsnow.o : LIS_historyMod.o
