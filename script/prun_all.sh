#!/bin/sh
if [ ! -f "../src/parallel" ]; then
    echo "compile parallel program..."
    g++ -std=c++11 -fopenmp ../src/parallel.cpp -o ../src/parallel
fi
for ((i=1; i<=12; i++)); do
    echo "case$i:-----------------------------------"
    ../src/parallel ../data/case$i.txt ../presult/res$i.txt
done