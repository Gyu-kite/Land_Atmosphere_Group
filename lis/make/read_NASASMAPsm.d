read_NASASMAPsm.o read_NASASMAPsm.d : read_NASASMAPsm.F90
read_NASASMAPsm.o : LIS_dataAssimMod.o
read_NASASMAPsm.o : LIS_mpiMod.o
read_NASASMAPsm.o : LIS_DAobservationsMod.o
read_NASASMAPsm.o : LIS_constantsMod.o
read_NASASMAPsm.o : LIS_pluginIndices.o
read_NASASMAPsm.o : NASASMAPsm_Mod.o
read_NASASMAPsm.o : LIS_logMod.o
read_NASASMAPsm.o : LIS_misc.h
read_NASASMAPsm.o : LIS_coreMod.o
read_NASASMAPsm.o : LIS_timeMgrMod.o
read_NASASMAPsm.o : map_utils.o
