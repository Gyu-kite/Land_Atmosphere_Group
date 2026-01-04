read_HYMAP_refq.o read_HYMAP_refq.d : read_HYMAP_refq.F90
read_HYMAP_refq.o : HYMAP_parmsMod.o
read_HYMAP_refq.o : LDT_logMod.o
read_HYMAP_refq.o : LDT_fileIOMod.o
read_HYMAP_refq.o : LDT_coreMod.o
