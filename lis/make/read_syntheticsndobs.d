read_syntheticsndobs.o read_syntheticsndobs.d : read_syntheticsndobs.F90
read_syntheticsndobs.o : LIS_DAobservationsMod.o
read_syntheticsndobs.o : LIS_misc.h
read_syntheticsndobs.o : LIS_coreMod.o
read_syntheticsndobs.o : LIS_logMod.o
read_syntheticsndobs.o : LIS_historyMod.o
read_syntheticsndobs.o : LIS_constantsMod.o
read_syntheticsndobs.o : LIS_NetCDF_inc.h
