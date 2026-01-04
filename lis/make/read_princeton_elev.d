read_princeton_elev.o read_princeton_elev.d : read_princeton_elev.F90
read_princeton_elev.o : LIS_metforcingMod.o
read_princeton_elev.o : LIS_logMod.o
read_princeton_elev.o : princeton_forcingMod.o
read_princeton_elev.o : LIS_coreMod.o
read_princeton_elev.o : LIS_fileIOMod.o
