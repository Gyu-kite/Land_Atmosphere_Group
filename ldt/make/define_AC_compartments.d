define_AC_compartments.o define_AC_compartments.d : define_AC_compartments.F90
define_AC_compartments.o : LDT_coreMod.o
define_AC_compartments.o : LDT_paramDataMod.o
define_AC_compartments.o : ac_utils.o
define_AC_compartments.o : AquaCrop_parmsMod.o
define_AC_compartments.o : LDT_constantsMod.o
define_AC_compartments.o : LDT_logMod.o
