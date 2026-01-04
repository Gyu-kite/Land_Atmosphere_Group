write_LPRM_AMSREsmobs.o write_LPRM_AMSREsmobs.d : write_LPRM_AMSREsmobs.F90
write_LPRM_AMSREsmobs.o : LIS_fileIOMod.o
write_LPRM_AMSREsmobs.o : LIS_logMod.o
write_LPRM_AMSREsmobs.o : LIS_coreMod.o
write_LPRM_AMSREsmobs.o : LIS_constantsMod.o
write_LPRM_AMSREsmobs.o : LIS_DAobservationsMod.o
write_LPRM_AMSREsmobs.o : LIS_historyMod.o
