#!/bin/sh
echo "compile serial program..."
g++ -std=c++11 -fopenmp ../src/serial.cpp -o ../src/serial
echo "compile parallel program..."
g++ -std=c++11 -fopenmp ../src/parallel.cpp -o ../src/parallel