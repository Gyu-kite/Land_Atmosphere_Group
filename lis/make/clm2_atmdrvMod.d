clm2_atmdrvMod.o clm2_atmdrvMod.d : clm2_atmdrvMod.F90
clm2_atmdrvMod.o : LIS_coreMod.o
clm2_atmdrvMod.o : clm2_varcon.o
clm2_atmdrvMod.o : LIS_precisionMod.o
clm2_atmdrvMod.o : LIS_metforcingMod.o
clm2_atmdrvMod.o : LIS_logMod.o
clm2_atmdrvMod.o : LIS_misc.h
clm2_atmdrvMod.o : LIS_FORC_AttributesMod.o
clm2_atmdrvMod.o : clm2_lsmMod.o
