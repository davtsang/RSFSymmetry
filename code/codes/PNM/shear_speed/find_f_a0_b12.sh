#!/bin/bash
for f in $(seq 4000000 250000 7000000)
do
for J in $(seq 25.0 1.4 36.2)
do
./find_f $f $J >"ctav_contours/J$J""_ctav$f.txt"
done
done
exit 0
