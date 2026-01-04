write_FLUXNETdata.o write_FLUXNETdata.d : write_FLUXNETdata.F90
write_FLUXNETdata.o : LIS_constantsMod.o
write_FLUXNETdata.o : LIS_coreMod.o
write_FLUXNETdata.o : LIS_historyMod.o
write_FLUXNETdata.o : LIS_fileIOMod.o
write_FLUXNETdata.o : LIS_logMod.o
