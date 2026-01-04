write_WUS_UCLAsnow.o write_WUS_UCLAsnow.d : write_WUS_UCLAsnow.F90
write_WUS_UCLAsnow.o : LIS_constantsMod.o
write_WUS_UCLAsnow.o : LIS_coreMod.o
write_WUS_UCLAsnow.o : LIS_DAobservationsMod.o
write_WUS_UCLAsnow.o : LIS_logMod.o
write_WUS_UCLAsnow.o : LIS_historyMod.o
write_WUS_UCLAsnow.o : LIS_fileIOMod.o
