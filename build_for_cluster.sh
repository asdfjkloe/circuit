#!/bin/zsh
g++ -std=c++14 -march=native -fpermissive -fopenmp -Ofast -fno-finite-math-only -I/home/fm093026/local_usr/include/ main.cpp -o circuit -L/home/fm093026/local_usr/lib -lopenblassingle -lgomp -lsuperlu -Wall
