climatology_forcingMod.o climatology_forcingMod.d : climatology_forcingMod.F90
climatology_forcingMod.o : LIS_timeMgrMod.o
climatology_forcingMod.o : LIS_constantsMod.o
climatology_forcingMod.o : climatology_SpatialInterpMod.o
climatology_forcingMod.o : climatology_VariablesMod.o
climatology_forcingMod.o : LIS_coreMod.o
climatology_forcingMod.o : LIS_misc.h
climatology_forcingMod.o : LIS_logMod.o
