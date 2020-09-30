#!/bin/bash
gcc -c -Wall -DNO_IEEE_INFINITY Time_2d.c
gfortran -c -Wall RIKsrf2.f90

gfortran -o RIKsrf2 -fopenmp RIKsrf2.o Time_2d.o 
