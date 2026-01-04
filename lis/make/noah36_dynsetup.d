noah36_dynsetup.o noah36_dynsetup.d : noah36_dynsetup.F90
noah36_dynsetup.o : module_sf_noah36lsm.o
noah36_dynsetup.o : LIS_coreMod.o
noah36_dynsetup.o : LIS_snowMod.o
noah36_dynsetup.o : noah36_lsmMod.o
noah36_dynsetup.o : LIS_timeMgrMod.o
noah36_dynsetup.o : LIS_albedoMod.o
noah36_dynsetup.o : LIS_emissMod.o
noah36_dynsetup.o : LIS_vegDataMod.o
