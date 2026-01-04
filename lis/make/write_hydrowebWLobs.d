write_hydrowebWLobs.o write_hydrowebWLobs.d : write_hydrowebWLobs.F90
write_hydrowebWLobs.o : LIS_historyMod.o
write_hydrowebWLobs.o : LIS_fileIOMod.o
write_hydrowebWLobs.o : LIS_coreMod.o
write_hydrowebWLobs.o : LIS_constantsMod.o
write_hydrowebWLobs.o : hydrowebWLobs_module.o
write_hydrowebWLobs.o : LIS_DAobservationsMod.o
write_hydrowebWLobs.o : LIS_logMod.o
