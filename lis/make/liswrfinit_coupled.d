liswrfinit_coupled.o liswrfinit_coupled.d : liswrfinit_coupled.F90
liswrfinit_coupled.o : LIS_surfaceModelMod.o
liswrfinit_coupled.o : LIS_metforcingMod.o
liswrfinit_coupled.o : LIS_paramsMod.o
liswrfinit_coupled.o : LIS_domainMod.o
liswrfinit_coupled.o : LIS_dataAssimMod.o
liswrfinit_coupled.o : LIS_irrigationMod.o
liswrfinit_coupled.o : LIS_DAobservationsMod.o
liswrfinit_coupled.o : LISWRFGridCompMod.o
liswrfinit_coupled.o : LIS_perturbMod.o
liswrfinit_coupled.o : LIS_tbotAdjustMod.o
liswrfinit_coupled.o : LIS_coreMod.o
liswrfinit_coupled.o : LIS_logMod.o
