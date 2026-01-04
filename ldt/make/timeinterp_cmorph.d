timeinterp_cmorph.o timeinterp_cmorph.d : timeinterp_cmorph.F90
timeinterp_cmorph.o : cmorph_forcingMod.o
timeinterp_cmorph.o : LDT_metforcingMod.o
timeinterp_cmorph.o : LDT_coreMod.o
timeinterp_cmorph.o : LDT_logMod.o
timeinterp_cmorph.o : LDT_FORC_AttributesMod.o
