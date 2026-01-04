write_PMW_snowobs.o write_PMW_snowobs.d : write_PMW_snowobs.F90
write_PMW_snowobs.o : LIS_logMod.o
write_PMW_snowobs.o : LIS_coreMod.o
write_PMW_snowobs.o : LIS_constantsMod.o
write_PMW_snowobs.o : LIS_fileIOMod.o
write_PMW_snowobs.o : LIS_historyMod.o
