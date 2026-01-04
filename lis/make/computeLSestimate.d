computeLSestimate.o computeLSestimate.d : computeLSestimate.F90
computeLSestimate.o : LSobjFunc_Mod.o
computeLSestimate.o : LIS_coreMod.o
computeLSestimate.o : LIS_optUEMod.o
computeLSestimate.o : LIS_mpiMod.o
computeLSestimate.o : LIS_logMod.o
