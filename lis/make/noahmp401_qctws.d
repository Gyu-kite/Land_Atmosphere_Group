noahmp401_qctws.o noahmp401_qctws.d : noahmp401_qctws.F90
noahmp401_qctws.o : NoahMP401_lsmMod.o
noahmp401_qctws.o : LIS_coreMod.o
noahmp401_qctws.o : LIS_logMod.o
noahmp401_qctws.o : module_sf_noahmplsm_401.o
