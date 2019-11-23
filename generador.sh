gfortran -g -c moduloArchivos.f90
gfortran -g moduloArchivos.o generadorDeCasos.f90 -o casos
./casos