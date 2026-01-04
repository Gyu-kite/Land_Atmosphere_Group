noahmp36_qc_soilmobs.o noahmp36_qc_soilmobs.d : noahmp36_qc_soilmobs.F90
noahmp36_qc_soilmobs.o : module_sf_noahlsm_36.o
noahmp36_qc_soilmobs.o : NoahMP36_lsmMod.o
noahmp36_qc_soilmobs.o : LIS_DAobservationsMod.o
noahmp36_qc_soilmobs.o : LIS_logMod.o
noahmp36_qc_soilmobs.o : LIS_constantsMod.o
noahmp36_qc_soilmobs.o : LIS_coreMod.o
