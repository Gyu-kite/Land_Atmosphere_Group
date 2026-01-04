read_AWRAL.o read_AWRAL.d : read_AWRAL.F90
read_AWRAL.o : LIS_misc.h
read_AWRAL.o : AWRAL_forcingMod.o
read_AWRAL.o : LIS_coreMod.o
read_AWRAL.o : LIS_constantsMod.o
read_AWRAL.o : LIS_logMod.o
read_AWRAL.o : LIS_metforcingMod.o
