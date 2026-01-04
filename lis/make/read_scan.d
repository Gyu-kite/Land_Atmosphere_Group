read_scan.o read_scan.d : read_scan.F90
read_scan.o : scan_forcingMod.o
read_scan.o : LIS_coreMod.o
read_scan.o : LIS_constantsMod.o
read_scan.o : LIS_metforcingMod.o
read_scan.o : LIS_logMod.o
