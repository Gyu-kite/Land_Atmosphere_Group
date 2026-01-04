read_CNRSdata_MPDI.o read_CNRSdata_MPDI.d : read_CNRSdata_MPDI.F90
read_CNRSdata_MPDI.o : LIS_logMod.o
read_CNRSdata_MPDI.o : LIS_mpiMod.o
read_CNRSdata_MPDI.o : LIS_fileIOMod.o
read_CNRSdata_MPDI.o : LIS_timeMgrMod.o
read_CNRSdata_MPDI.o : LIS_coreMod.o
read_CNRSdata_MPDI.o : map_utils.o
read_CNRSdata_MPDI.o : LIS_RTMMod.o
read_CNRSdata_MPDI.o : CNRS_em_obsMod.o
