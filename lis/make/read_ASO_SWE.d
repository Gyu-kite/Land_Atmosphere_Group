read_ASO_SWE.o read_ASO_SWE.d : read_ASO_SWE.F90
read_ASO_SWE.o : LIS_DAobservationsMod.o
read_ASO_SWE.o : LIS_coreMod.o
read_ASO_SWE.o : ASO_SWE_Mod.o
read_ASO_SWE.o : map_utils.o
read_ASO_SWE.o : LIS_timeMgrMod.o
read_ASO_SWE.o : LIS_constantsMod.o
read_ASO_SWE.o : LIS_misc.h
read_ASO_SWE.o : LIS_logMod.o
read_ASO_SWE.o : LIS_mpiMod.o
read_ASO_SWE.o : UTM_utils.o
read_ASO_SWE.o : LIS_dataAssimMod.o
