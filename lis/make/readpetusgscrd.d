readpetusgscrd.o readpetusgscrd.d : readpetusgscrd.F90
readpetusgscrd.o : petusgs_forcingMod.o
readpetusgscrd.o : LIS_logMod.o
readpetusgscrd.o : LIS_coreMod.o
