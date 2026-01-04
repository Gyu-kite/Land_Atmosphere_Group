read_ecmwf_elev.o read_ecmwf_elev.d : read_ecmwf_elev.F90
read_ecmwf_elev.o : LIS_metforcingMod.o
read_ecmwf_elev.o : LIS_fileIOMod.o
read_ecmwf_elev.o : LIS_logMod.o
read_ecmwf_elev.o : ecmwf_forcingMod.o
read_ecmwf_elev.o : LIS_constantsMod.o
read_ecmwf_elev.o : LIS_coreMod.o
