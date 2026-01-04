ac72_main.o ac72_main.d : ac72_main.F90
ac72_main.o : global.o
ac72_main.o : run.o
ac72_main.o : ac72_prep_f.o
ac72_main.o : LIS_histDataMod.o
ac72_main.o : LIS_constantsMod.o
ac72_main.o : LIS_logMod.o
ac72_main.o : startunit.o
ac72_main.o : kinds.o
ac72_main.o : ac72_lsmMod.o
ac72_main.o : project_input.o
ac72_main.o : LIS_coreMod.o
ac72_main.o : LIS_timeMgrMod.o
