#! /bin/bash
hilos=$(nproc)
gfortran -g -c -fcoarray=lib RBGS.f90 -lcaf_mpi 
gfortran -g -c SELs.f90
gfortran -g -c moduloArchivos.f90
gfortran -g -c arreglos.f90
gfortran -g moduloArchivos.o arreglos.o SELs.o RBGS.o -fcoarray=lib $1 -lcaf_mpi -o $(basename $1 .f90)
cafrun -np $hilos ./$(basename $1 .f90)