write_ARMdata.o write_ARMdata.d : write_ARMdata.F90
write_ARMdata.o : LIS_coreMod.o
write_ARMdata.o : LIS_historyMod.o
write_ARMdata.o : LIS_logMod.o
write_ARMdata.o : LIS_fileIOMod.o
write_ARMdata.o : LIS_constantsMod.o
