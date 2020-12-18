#!/bin/bash
for K in $(seq 21 -1 00)
do
for L in $(seq 20.01 10.00 40.01)
do
for J in $(seq 26.0 1.0 35.0)
do
./TOV $J $L $K 1.40
done
done
done

for K in $(seq 21 -1 00)
do
for L in $(seq 50.02 10.00 70.02)
do
for J in $(seq 26.0 1.0 35.0)
do
./TOV $J $L $K 1.40
done
done
done

for K in $(seq 21 -1 00)
do
for L in $(seq 80.03 10.00 100.03)
do
for J in $(seq 26.0 1.0 35.0)
do
./TOV $J $L $K 1.40
done
done
done
exit 0





