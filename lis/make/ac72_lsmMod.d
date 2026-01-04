ac72_lsmMod.o ac72_lsmMod.d : ac72_lsmMod.F90
ac72_lsmMod.o : ac72_module.o
ac72_lsmMod.o : LIS_timeMgrMod.o
ac72_lsmMod.o : LIS_logMod.o
ac72_lsmMod.o : LIS_constantsMod.o
ac72_lsmMod.o : LIS_coreMod.o
ac72_lsmMod.o : LIS_surfaceModelDataMod.o
