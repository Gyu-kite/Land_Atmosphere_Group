read_WUS_UCLAsnow.o read_WUS_UCLAsnow.d : read_WUS_UCLAsnow.F90
read_WUS_UCLAsnow.o : WUS_UCLAsnowMod.o
read_WUS_UCLAsnow.o : LIS_pluginIndices.o
read_WUS_UCLAsnow.o : map_utils.o
read_WUS_UCLAsnow.o : LIS_dataAssimMod.o
read_WUS_UCLAsnow.o : LIS_mpiMod.o
read_WUS_UCLAsnow.o : LIS_constantsMod.o
read_WUS_UCLAsnow.o : LIS_timeMgrMod.o
read_WUS_UCLAsnow.o : LIS_misc.h
read_WUS_UCLAsnow.o : LIS_coreMod.o
read_WUS_UCLAsnow.o : LIS_DAobservationsMod.o
read_WUS_UCLAsnow.o : LIS_logMod.o
