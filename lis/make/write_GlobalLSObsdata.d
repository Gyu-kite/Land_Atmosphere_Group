write_GlobalLSObsdata.o write_GlobalLSObsdata.d : write_GlobalLSObsdata.F90
write_GlobalLSObsdata.o : LIS_constantsMod.o
write_GlobalLSObsdata.o : LIS_logMod.o
write_GlobalLSObsdata.o : LIS_fileIOMod.o
write_GlobalLSObsdata.o : LIS_historyMod.o
write_GlobalLSObsdata.o : LIS_coreMod.o
