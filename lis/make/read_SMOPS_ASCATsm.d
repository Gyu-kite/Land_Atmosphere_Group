read_SMOPS_ASCATsm.o read_SMOPS_ASCATsm.d : read_SMOPS_ASCATsm.F90
read_SMOPS_ASCATsm.o : LIS_mpiMod.o
read_SMOPS_ASCATsm.o : LIS_dataAssimMod.o
read_SMOPS_ASCATsm.o : LIS_coreMod.o
read_SMOPS_ASCATsm.o : LIS_timeMgrMod.o
read_SMOPS_ASCATsm.o : map_utils.o
read_SMOPS_ASCATsm.o : LIS_misc.h
read_SMOPS_ASCATsm.o : LIS_constantsMod.o
read_SMOPS_ASCATsm.o : SMOPS_ASCATsm_Mod.o
read_SMOPS_ASCATsm.o : LIS_pluginIndices.o
read_SMOPS_ASCATsm.o : LIS_DAobservationsMod.o
read_SMOPS_ASCATsm.o : LIS_logMod.o
