read_nam242_elev.o read_nam242_elev.d : read_nam242_elev.F90
read_nam242_elev.o : LIS_logMod.o
read_nam242_elev.o : LIS_coreMod.o
read_nam242_elev.o : nam242_forcingMod.o
read_nam242_elev.o : LIS_metforcingMod.o
read_nam242_elev.o : LIS_fileIOMod.o
