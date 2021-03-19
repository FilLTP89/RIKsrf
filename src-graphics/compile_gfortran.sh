#!/bin/bash
gfortran -c -Wall RIKsrf2anime.f90
gfortran -o RIKsrf2anime -fopenmp RIKsrf2anime.o
