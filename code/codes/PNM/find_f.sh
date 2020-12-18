#!/bin/bash
for f in $(seq 100 20 300)
do
for J in $(seq 25.0 1.4 36.2)
do
./find_f $f $J 0 0 >"contours_a0_b12/f$f""_J$J.txt"
done
done
exit 0
