read_petusgs.o read_petusgs.d : read_petusgs.F90
read_petusgs.o : LIS_metforcingMod.o
read_petusgs.o : LIS_logMod.o
read_petusgs.o : fbil_module.o
read_petusgs.o : petusgs_forcingMod.o
read_petusgs.o : LIS_coreMod.o
