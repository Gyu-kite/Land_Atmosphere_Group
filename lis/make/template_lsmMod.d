template_lsmMod.o template_lsmMod.d : template_lsmMod.F90
template_lsmMod.o : LIS_misc.h
template_lsmMod.o : LIS_surfaceModelDataMod.o
template_lsmMod.o : LIS_coreMod.o
template_lsmMod.o : template_module.o
template_lsmMod.o : LIS_logMod.o
template_lsmMod.o : LIS_timeMgrMod.o
