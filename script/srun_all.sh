#!/bin/sh
if [ ! -f "../src/serial" ]; then
    echo "compile serial program..."
    g++ -std=c++11 -fopenmp ../src/serial.cpp -o ../src/serial
fi
for ((i=1; i<=12; i++)); do
    echo "case$i:-----------------------------------"
    ../src/serial ../data/case$i.txt ../sresult/res$i.txt
done