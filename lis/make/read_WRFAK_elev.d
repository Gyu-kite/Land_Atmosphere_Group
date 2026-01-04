read_WRFAK_elev.o read_WRFAK_elev.d : read_WRFAK_elev.F90
read_WRFAK_elev.o : LIS_fileIOMod.o
read_WRFAK_elev.o : LIS_coreMod.o
read_WRFAK_elev.o : LIS_logMod.o
read_WRFAK_elev.o : LIS_metforcingMod.o
