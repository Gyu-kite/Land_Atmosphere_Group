NLDAS_routing_run.o NLDAS_routing_run.d : NLDAS_routing_run.F90
NLDAS_routing_run.o : LIS_coreMod.o
NLDAS_routing_run.o : NLDAS_routingMod.o
NLDAS_routing_run.o : LIS_routingMod.o
NLDAS_routing_run.o : LIS_historyMod.o
NLDAS_routing_run.o : LIS_histDataMod.o
NLDAS_routing_run.o : LIS_logMod.o
NLDAS_routing_run.o : LIS_timeMgrMod.o
