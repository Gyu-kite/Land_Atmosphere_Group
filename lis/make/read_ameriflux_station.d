read_ameriflux_station.o read_ameriflux_station.d : read_ameriflux_station.F90
read_ameriflux_station.o : LIS_logMod.o
read_ameriflux_station.o : AmerifluxobsMod.o
read_ameriflux_station.o : LIS_coreMod.o
read_ameriflux_station.o : LIS_constantsMod.o
read_ameriflux_station.o : LIS_timeMgrMod.o
