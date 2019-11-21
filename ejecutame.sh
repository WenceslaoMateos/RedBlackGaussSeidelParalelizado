#! /bin/bash
hilos=$(nproc)
gfortran -g -c -fcoarray=lib RBGS.f90 -lcaf_mpi 
gfortran -g -c SELs.f90
gfortran -g -c arreglos.f90
gfortran -g arreglos.o SELs.o RBGS.o -fcoarray=lib main.f90 -lcaf_mpi -o main
cafrun -np $hilos ./main