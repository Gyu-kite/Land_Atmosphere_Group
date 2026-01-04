write_syntheticsndobs.o write_syntheticsndobs.d : write_syntheticsndobs.F90
write_syntheticsndobs.o : LIS_fileIOMod.o
write_syntheticsndobs.o : LIS_logMod.o
write_syntheticsndobs.o : LIS_historyMod.o
write_syntheticsndobs.o : LIS_constantsMod.o
write_syntheticsndobs.o : LIS_coreMod.o
write_syntheticsndobs.o : LIS_DAobservationsMod.o
