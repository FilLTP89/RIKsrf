#!/bin/bash
gcc -c -Wall -DNO_IEEE_INFINITY Time_2d.c
mpif90 -c -Wall -g -fbacktrace RIKsrf_parallel.f90

mpif90 -o RIKsrf_parallel -fopenmp RIKsrf2.o Time_2d.o 
