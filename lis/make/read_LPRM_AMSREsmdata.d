read_LPRM_AMSREsmdata.o read_LPRM_AMSREsmdata.d : read_LPRM_AMSREsmdata.F90
read_LPRM_AMSREsmdata.o : LIS_coreMod.o
read_LPRM_AMSREsmdata.o : LIS_fileIOMod.o
read_LPRM_AMSREsmdata.o : LPRM_AMSREsm_obsMod.o
read_LPRM_AMSREsmdata.o : LIS_misc.h
read_LPRM_AMSREsmdata.o : map_utils.o
read_LPRM_AMSREsmdata.o : LIS_logMod.o
read_LPRM_AMSREsmdata.o : LIS_constantsMod.o
read_LPRM_AMSREsmdata.o : LIS_mpiMod.o
read_LPRM_AMSREsmdata.o : LIS_timeMgrMod.o
