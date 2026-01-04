read_gdaselev_ldtproc.o read_gdaselev_ldtproc.d : read_gdaselev_ldtproc.F90
read_gdaselev_ldtproc.o : gdas_forcingMod.o
read_gdaselev_ldtproc.o : LDT_logMod.o
read_gdaselev_ldtproc.o : LDT_fileIOMod.o
read_gdaselev_ldtproc.o : LDT_coreMod.o
read_gdaselev_ldtproc.o : LDT_metforcingMod.o
