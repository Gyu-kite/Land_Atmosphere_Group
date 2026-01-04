run.o run.d : run.f90
run.o : preparefertilitysalinity.o
run.o : utils.o
run.o : simul.o
run.o : project_input.o
run.o : inforesults.o
run.o : tempprocessing.o
run.o : rootunit.o
run.o : kinds.o
run.o : global.o
run.o : climprocessing.o
